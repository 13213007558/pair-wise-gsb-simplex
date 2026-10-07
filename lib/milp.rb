# Branch-and-bound solver for mixed-integer linear programs.
#
# Builds on the same problem format as Simplex (maximization, <= constraints,
# nonnegative variables) and adds per-variable integrality:
#
#   milp = Milp.new(
#     [3, 2],            # objective coefficients
#     [[2, 1], [1, 2]],  # constraint matrix
#     [4, 3],            # right hand side
#     [:integer, :binary]
#   )
#   result = milp.solve
#   result.status    # => :optimal
#   result.solution  # best integer feasible solution found
#   result.objective # its objective value
#   result.bound     # global bound on the optimal objective
#   result.nodes     # number of LP relaxations solved
#
# The LP relaxation at each node is solved with an exact (Rational arithmetic)
# two-phase simplex, so results are never produced by rounding the continuous
# optimum.
class Milp
  DEFAULT_MAX_NODES = 10_000

  VARIABLE_TYPES    = [:continuous, :integer, :binary].freeze
  BRANCH_STRATEGIES = [:first, :most_fractional].freeze
  NODE_SELECTIONS   = [:depth_first, :best_bound].freeze

  INFINITY = Float::INFINITY

  # status:   :optimal              - tree exhausted, +solution+ is proven optimal
  #           :no_integer_solution  - relaxation feasible, but no integer feasible point exists
  #           :infeasible           - the LP relaxation itself is infeasible
  #           :unbounded            - the LP relaxation is unbounded
  #           :node_limit_exceeded  - node budget exhausted before optimality was proven
  # solution: best integer feasible point found, or nil
  # objective: objective value of +solution+, or nil
  # bound:    global upper bound on the optimal objective value, or nil when unknown
  # nodes:    number of LP relaxations solved
  Result = Struct.new(:status, :solution, :objective, :bound, :nodes) do
    def optimal?
      status == :optimal
    end
  end

  Node = Struct.new(:lower, :upper, :bound, :sequence)

  # c, a, b:        same format as Simplex (max c.x s.t. a.x <= b, x >= 0)
  # variable_types: one of :continuous, :integer, :binary per variable
  #                 (defaults to all :continuous, which reduces to plain LP)
  # max_nodes:      budget of LP relaxations to solve
  # branch_strategy: :first (fixed variable order) or :most_fractional
  # node_selection:  :depth_first or :best_bound
  def initialize(c, a, b, variable_types = nil,
                 max_nodes: DEFAULT_MAX_NODES,
                 branch_strategy: :most_fractional,
                 node_selection: :best_bound)
    num_vars = c.length
    unless a.all? { |row| row.length == num_vars } && b.length == a.length
      raise ArgumentError, "Input arrays have mismatched dimensions"
    end
    variable_types ||= [:continuous] * num_vars
    unless variable_types.length == num_vars &&
           (variable_types - VARIABLE_TYPES).empty?
      raise ArgumentError, "variable_types must contain one of " \
                           "#{VARIABLE_TYPES.inspect} per variable"
    end
    unless BRANCH_STRATEGIES.include?(branch_strategy)
      raise ArgumentError, "unknown branch strategy #{branch_strategy.inspect}"
    end
    unless NODE_SELECTIONS.include?(node_selection)
      raise ArgumentError, "unknown node selection #{node_selection.inspect}"
    end
    unless max_nodes.is_a?(Integer) && max_nodes >= 0
      raise ArgumentError, "max_nodes must be a non-negative integer"
    end

    @c = c
    @a = a
    @b = b
    @num_vars = num_vars
    @variable_types = variable_types
    @integer_vars = (0...num_vars).select { |j| variable_types[j] != :continuous }
    @max_nodes = max_nodes
    @branch_strategy = branch_strategy
    @node_selection = node_selection
  end

  def solve
    @nodes = 0
    @incumbent = nil
    @incumbent_objective = nil
    @sequence = 0
    root_infeasible = false

    open_nodes = [make_node([0] * @num_vars, initial_uppers, INFINITY)]

    until open_nodes.empty?
      return node_limit_result(open_nodes) if @nodes >= @max_nodes

      node = select_node(open_nodes)
      relaxation = solve_relaxation(node)

      case relaxation.status
      when :infeasible
        root_infeasible = true if @nodes == 1
        next
      when :unbounded
        # Only possible at the root: child nodes are more constrained.
        return Result.new(:unbounded, nil, nil, nil, @nodes)
      end

      next if @incumbent_objective && relaxation.objective <= @incumbent_objective

      if integer_solution?(relaxation.x)
        @incumbent = normalize(relaxation.x)
        @incumbent_objective = normalize_value(relaxation.objective)
      else
        branch(node, relaxation).reverse_each { |child| open_nodes << child }
      end
    end

    if @incumbent
      Result.new(:optimal, @incumbent, @incumbent_objective,
                 @incumbent_objective, @nodes)
    elsif root_infeasible
      Result.new(:infeasible, nil, nil, nil, @nodes)
    else
      Result.new(:no_integer_solution, nil, nil, nil, @nodes)
    end
  end

  private

  def make_node(lower, upper, bound)
    @sequence += 1
    Node.new(lower, upper, bound, @sequence)
  end

  def initial_uppers
    @variable_types.map { |type| type == :binary ? 1 : INFINITY }
  end

  def solve_relaxation(node)
    @nodes += 1
    a = @a.map(&:dup)
    b = @b.dup
    @num_vars.times do |j|
      if node.upper[j] < INFINITY
        row = Array.new(@num_vars) { 0 }
        row[j] = 1
        a << row
        b << node.upper[j]
      end
      if node.lower[j] > 0
        row = Array.new(@num_vars) { 0 }
        row[j] = -1
        a << row
        b << -node.lower[j]
      end
    end
    LP.solve(@c, a, b)
  end

  def select_node(open_nodes)
    case @node_selection
    when :depth_first
      open_nodes.pop
    when :best_bound
      best = open_nodes.max_by { |node| [node.bound, node.sequence] }
      open_nodes.delete(best)
      best
    end
  end

  def integer_solution?(x)
    @integer_vars.all? { |j| integral?(x[j]) }
  end

  def integral?(value)
    value == value.floor
  end

  # Deterministic: :first picks the lowest-index fractional integer variable,
  # :most_fractional picks the one closest to 1/2 (ties broken by index).
  def branch_variable(x)
    fractional = @integer_vars.reject { |j| integral?(x[j]) }
    case @branch_strategy
    when :first
      fractional.first
    when :most_fractional
      fractional.max_by do |j|
        part = x[j] - x[j].floor
        [[part, 1 - part].min, -j]
      end
    end
  end

  # Children are created down-branch (x_j <= floor) first, then up-branch
  # (x_j >= ceil); callers push them so that the down branch is explored
  # first under :depth_first.
  def branch(node, relaxation)
    j = branch_variable(relaxation.x)
    value = relaxation.x[j]

    down = make_node(node.lower.dup, node.upper.dup, relaxation.objective)
    down.upper[j] = value.floor

    up = make_node(node.lower.dup, node.upper.dup, relaxation.objective)
    up.lower[j] = value.floor + 1

    [down, up]
  end

  def node_limit_result(open_nodes)
    bounds = open_nodes.map(&:bound).reject { |bound| bound == INFINITY }
    bounds << @incumbent_objective if @incumbent_objective
    Result.new(:node_limit_exceeded, @incumbent, @incumbent_objective,
               bounds.max, @nodes)
  end

  def normalize(x)
    x.map { |value| normalize_value(value) }
  end

  def normalize_value(value)
    integral?(value) ? value.floor : value
  end

  # Exact two-phase primal simplex for  max c.x  s.t.  A.x <= b, x >= 0
  # with arbitrary signs in b. All arithmetic is Rational, so feasibility
  # and integrality decisions are exact. Uses Bland's rule (no cycling).
  class LP
    Result = Struct.new(:status, :x, :objective)

    def self.solve(c, a, b)
      new(c, a, b).solve
    end

    def initialize(c, a, b)
      @num_vars = c.length
      @num_constraints = b.length
      # Columns: original variables, slack variables, one artificial variable.
      @num_columns = @num_vars + @num_constraints + 1
      @artificial = @num_columns - 1

      @tableau = Array.new(@num_constraints) do |i|
        row = Array.new(@num_columns + 1) { Rational(0) }
        @num_vars.times { |j| row[j] = a[i][j].to_r }
        row[@num_vars + i] = Rational(1)
        row[@artificial] = Rational(-1)
        row[-1] = b[i].to_r
        row
      end
      @basis = (0...@num_constraints).map { |i| @num_vars + i }
      @objective = Array.new(@num_columns + 1) { Rational(0) }
      @original_c = c
    end

    def solve
      # Phase 1: maximize -(artificial). If some b_i < 0 the slack basis is
      # infeasible, so first pivot the artificial variable into the row with
      # the most negative rhs, which restores primal feasibility.
      set_objective(phase_one_costs)
      most_negative = (0...@num_constraints).min_by { |i| @tableau[i][-1] }
      if most_negative && @tableau[most_negative][-1] < 0
        pivot(most_negative, @artificial)
      end
      run_simplex

      return Result.new(:infeasible, nil, nil) if objective_value < 0

      drive_out_artificial
      drop_artificial_column

      # Phase 2: optimize the original objective from the feasible basis.
      set_objective(original_costs)
      return Result.new(:unbounded, nil, nil) if run_simplex == :unbounded

      x = Array.new(@num_vars) { Rational(0) }
      @basis.each_with_index do |var, i|
        x[var] = @tableau[i][-1] if var < @num_vars
      end
      objective = Rational(0)
      @num_vars.times { |j| objective += @original_c[j].to_r * x[j] }
      Result.new(:optimal, x, objective)
    end

    private

    def phase_one_costs
      costs = Array.new(@num_columns) { Rational(0) }
      costs[@artificial] = Rational(-1)
      costs
    end

    def original_costs
      (0...@num_columns).map do |j|
        j < @num_vars ? @original_c[j].to_r : Rational(0)
      end
    end

    # Rebuilds the objective row (reduced costs, last entry = -z) for the
    # current basis from the given cost vector.
    def set_objective(costs)
      @num_columns.times do |j|
        reduced = costs[j]
        @num_constraints.times do |i|
          reduced -= costs[@basis[i]] * @tableau[i][j]
        end
        @objective[j] = reduced
      end
      value = Rational(0)
      @num_constraints.times do |i|
        value += costs[@basis[i]] * @tableau[i][-1]
      end
      @objective[-1] = -value
    end

    def objective_value
      -@objective[-1]
    end

    def run_simplex
      loop do
        entering = (0...@num_columns).find { |j| @objective[j] > 0 }
        return :optimal unless entering

        pivot_row = nil
        best_ratio = nil
        @num_constraints.times do |i|
          coeff = @tableau[i][entering]
          next unless coeff > 0
          ratio = @tableau[i][-1] / coeff
          if best_ratio.nil? || ratio < best_ratio ||
             (ratio == best_ratio && @basis[i] < @basis[pivot_row])
            best_ratio = ratio
            pivot_row = i
          end
        end
        return :unbounded unless pivot_row

        pivot(pivot_row, entering)
      end
    end

    def pivot(row_ix, column_ix)
      pivot_row = @tableau[row_ix]
      ratio = Rational(1, pivot_row[column_ix])
      pivot_row.size.times { |j| pivot_row[j] *= ratio }

      @num_constraints.times do |i|
        next if i == row_ix
        eliminate(@tableau[i], pivot_row, column_ix)
      end
      eliminate(@objective, pivot_row, column_ix)

      @basis[row_ix] = column_ix
    end

    def eliminate(row, pivot_row, column_ix)
      factor = row[column_ix]
      return if factor.zero?
      row.size.times { |j| row[j] -= factor * pivot_row[j] }
    end

    # If the artificial variable is still basic its value is 0, so pivot it
    # out on any nonzero non-artificial entry; if its row is all zeros the
    # constraint is redundant and can be dropped.
    def drive_out_artificial
      @num_constraints.times do |i|
        next unless @basis[i] == @artificial
        column = (0...@artificial).find { |j| !@tableau[i][j].zero? }
        if column
          pivot(i, column)
        else
          @tableau.delete_at(i)
          @basis.delete_at(i)
          @num_constraints -= 1
          redo
        end
      end
    end

    def drop_artificial_column
      @tableau.each { |row| row.delete_at(@artificial) }
      @objective.delete_at(@artificial)
      @num_columns -= 1
    end
  end
end

