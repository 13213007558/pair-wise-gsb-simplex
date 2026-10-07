require 'simplex'

class Simplex
  # Branch and bound for (mixed-)integer linear programs, built on LP
  # relaxations. The continuous optimum is never rounded: an integer
  # solution is only recorded when a relaxation optimum already satisfies
  # every integrality constraint exactly (all arithmetic is rational).
  module MILP
    VARIABLE_TYPES    = [:continuous, :integer, :binary].freeze
    NODE_STRATEGIES   = [:best_bound, :depth_first].freeze
    BRANCH_STRATEGIES = [:most_fractional, :first].freeze
    DEFAULT_MAX_NODES = 10_000

    # Exact two-phase simplex for max c^T x subject to A x <= b, x >= 0.
    # Unlike Simplex it accepts right-hand sides of any sign, which branch
    # constraints such as x_j >= ceil(v) require.
    class RelaxedLP
      DEFAULT_MAX_PIVOTS = 100_000

      attr_reader :status, :solution, :objective_value

      def initialize(c, a, b, max_pivots: DEFAULT_MAX_PIVOTS)
        @num_vars    = c.length
        @c           = c.map {|v| Rational(v) }
        @a           = a.map {|row| row.map {|v| Rational(v) } }
        @b           = b.map {|v| Rational(v) }
        @max_pivots  = max_pivots
        @pivot_count = 0
      end

      def solve
        artificial = @num_vars
        total      = @num_vars + 1 + @b.length

        rows = (0...@b.length).map do |i|
          row = Array.new(total, 0)
          @num_vars.times {|j| row[j] = @a[i][j] }
          row[artificial] = -1
          row[@num_vars + 1 + i] = 1
          row
        end
        rhs   = @b.dup
        basis = (0...@b.length).map {|i| @num_vars + 1 + i }

        # Phase 1: maximize -artificial to find a feasible basis.
        costs = Array.new(total, 0)
        costs[artificial] = -1
        z = 0
        unless rhs.empty?
          most_violated = rhs.each_index.min_by {|i| rhs[i] }
          if rhs[most_violated] < 0
            z = pivot(rows, rhs, basis, costs, z, most_violated, artificial)
          end
        end
        z, = optimize(rows, rhs, basis, costs, z, artificial)
        if z < 0
          @status = :infeasible
          return self
        end

        # Drive the artificial variable out of the basis (it is zero).
        if (row_ix = basis.index(artificial))
          col = (0...total).find do |j|
            j != artificial && !basis.include?(j) && rows[row_ix][j] != 0
          end
          if col
            pivot(rows, rhs, basis, costs, z, row_ix, col)
          else
            rows.delete_at(row_ix)
            rhs.delete_at(row_ix)
            basis.delete_at(row_ix)
          end
        end

        # Phase 2: restore the real objective and re-optimize.
        objective = @c + Array.new(total - @num_vars, 0)
        costs = objective.dup
        z = 0
        basis.each_with_index do |basic, i|
          next if objective[basic] == 0
          z += objective[basic] * rhs[i]
          total.times {|j| costs[j] -= objective[basic] * rows[i][j] }
        end
        z, status = optimize(rows, rhs, basis, costs, z, artificial)
        if status == :unbounded
          @status = :unbounded
          return self
        end

        @status = :optimal
        @objective_value = z
        @solution = Array.new(@num_vars, 0)
        basis.each_with_index do |basic, i|
          @solution[basic] = rhs[i] if basic < @num_vars
        end
        self
      end

      private

      # Bland's rule (lowest-index entering and leaving variable) keeps the
      # exact rational arithmetic free of cycling.
      def optimize(rows, rhs, basis, costs, z, excluded)
        loop do
          entering = nil
          costs.each_with_index do |cost, j|
            next if j == excluded || basis.include?(j)
            if cost > 0
              entering = j
              break
            end
          end
          return [z, :optimal] unless entering

          leaving    = nil
          best_ratio = nil
          rows.each_with_index do |row, i|
            coeff = row[entering]
            next unless coeff > 0
            ratio = rhs[i] / coeff
            if best_ratio.nil? || ratio < best_ratio ||
               (ratio == best_ratio && basis[i] < basis[leaving])
              best_ratio = ratio
              leaving    = i
            end
          end
          return [z, :unbounded] unless leaving

          @pivot_count += 1
          raise "Too many pivots" if @pivot_count > @max_pivots
          z = pivot(rows, rhs, basis, costs, z, leaving, entering)
        end
      end

      def pivot(rows, rhs, basis, costs, z, row_ix, col_ix)
        pivot_row   = rows[row_ix]
        pivot_value = pivot_row[col_ix]
        if pivot_value != 1
          pivot_row.map! {|v| v / pivot_value }
          rhs[row_ix] /= pivot_value
        end
        rows.each_with_index do |row, i|
          next if i == row_ix
          coeff = row[col_ix]
          next if coeff == 0
          pivot_row.each_with_index {|v, j| row[j] -= coeff * v }
          rhs[i] -= coeff * rhs[row_ix]
        end
        cost = costs[col_ix]
        if cost != 0
          pivot_row.each_with_index {|v, j| costs[j] -= cost * v }
          z += cost * rhs[row_ix]
        end
        basis[row_ix] = col_ix
        z
      end
    end

    # Outcome of a branch-and-bound run.
    #
    # status is one of:
    #   :optimal             - search tree exhausted with an incumbent; the
    #                          solution is proven optimal
    #   :infeasible          - the LP relaxation itself is infeasible
    #   :no_integer_solution - continuous relaxation feasible, but the
    #                          exhausted search found no integer solution
    #   :unbounded           - the LP relaxation is unbounded
    #   :node_limit          - the node budget ran out; any incumbent is
    #                          returned together with a valid global bound
    class Result
      attr_reader :status, :solution, :objective, :bound, :nodes

      def initialize(status:, solution:, objective:, bound:, nodes:)
        @status    = status
        @solution  = solution
        @objective = objective
        @bound     = bound
        @nodes     = nodes
      end

      def optimal?
        @status == :optimal
      end
    end

    # Maximizing branch and bound. Deterministic: node and branching
    # choices are fully decided by the configured strategies, with
    # index/sequence tie-breaks.
    class BranchAndBound
      Node = Struct.new(:lower, :upper, :bound, :seq)

      def initialize(c, a, b, types, max_nodes:, node_strategy:, branch_strategy:)
        @c               = c.map {|v| Rational(v) }
        @a               = a.map {|row| row.map {|v| Rational(v) } }
        @b               = b.map {|v| Rational(v) }
        @types           = types
        @max_nodes       = max_nodes
        @node_strategy   = node_strategy
        @branch_strategy = branch_strategy
        @num_vars        = @c.length
        @integer_vars    = (0...@num_vars).select {|j| @types[j] != :continuous }
      end

      def solve
        nodes         = 0
        incumbent     = nil
        incumbent_obj = nil
        open_nodes    = [Node.new({}, {}, nil, 0)]
        next_seq      = 1
        status        = nil

        until open_nodes.empty?
          if nodes >= @max_nodes
            status = :node_limit
            break
          end

          node = select_node(open_nodes)
          lp   = RelaxedLP.new(@c, *node_constraints(node)).solve
          nodes += 1

          if lp.status == :infeasible
            status = :infeasible if nodes == 1
            next
          end
          if lp.status == :unbounded
            # Can only happen at the root: children restrict the region.
            return Result.new(status: :unbounded, solution: nil,
                              objective: nil, bound: nil, nodes: nodes)
          end

          obj = lp.objective_value
          next if incumbent_obj && obj <= incumbent_obj

          x = lp.solution
          if integer_feasible?(x)
            incumbent     = x
            incumbent_obj = obj
            open_nodes.reject! {|n| n.bound <= incumbent_obj }
            next
          end

          var   = branch_variable(x)
          value = x[var]
          # The down branch is enqueued last so that it is explored first.
          up   = Node.new(node.lower.merge(var => value.ceil), node.upper.dup,
                          obj, next_seq)
          down = Node.new(node.lower.dup, node.upper.merge(var => value.floor),
                          obj, next_seq + 1)
          next_seq += 2
          open_nodes.push(up, down)
        end

        status ||= incumbent ? :optimal : :no_integer_solution
        bound = case status
                when :optimal
                  incumbent_obj
                when :node_limit
                  (open_nodes.map(&:bound) + [incumbent_obj]).compact.max
                end
        Result.new(status: status,
                   solution: incumbent && incumbent.map {|v| normalize(v) },
                   objective: incumbent_obj && normalize(incumbent_obj),
                   bound: bound && normalize(bound),
                   nodes: nodes)
      end

      private

      def select_node(open_nodes)
        node =
          case @node_strategy
          when :best_bound
            # Largest relaxation bound first; ties go to the most recently
            # created node, so the down branch is explored first.
            open_nodes.max_by {|n| [n.bound, n.seq] }
          when :depth_first
            open_nodes.last
          end
        open_nodes.delete(node)
      end

      def node_constraints(node)
        a = @a.dup
        b = @b.dup
        @types.each_with_index do |type, j|
          next unless type == :binary
          row = Array.new(@num_vars, 0)
          row[j] = 1
          a << row
          b << 1
        end
        node.upper.each do |j, upper|
          row = Array.new(@num_vars, 0)
          row[j] = 1
          a << row
          b << upper
        end
        node.lower.each do |j, lower|
          row = Array.new(@num_vars, 0)
          row[j] = -1
          a << row
          b << -lower
        end
        [a, b]
      end

      def integer_feasible?(x)
        @integer_vars.all? {|j| x[j].denominator == 1 }
      end

      def branch_variable(x)
        fractional = @integer_vars.select {|j| x[j].denominator != 1 }
        case @branch_strategy
        when :first
          fractional.first
        when :most_fractional
          fractional.max_by do |j|
            frac = x[j] - x[j].floor
            [[frac, 1 - frac].min, -j]
          end
        end
      end

      def normalize(value)
        value.denominator == 1 ? value.to_i : value
      end
    end
  end

  # Solves the problem as a (mixed-)integer program by branch and bound on
  # LP relaxations.
  #
  # variable_types: one of :continuous, :integer or :binary per variable
  #   (may also be passed to Simplex.new as variable_types:).
  # max_nodes: node budget; the run stops with status :node_limit once this
  #   many LP relaxations have been solved.
  # node_strategy: :best_bound (default) or :depth_first.
  # branch_strategy: :most_fractional (default) or :first.
  #
  # Returns a MILP::Result; see there for the possible statuses.
  def solve_milp(variable_types: nil, max_nodes: MILP::DEFAULT_MAX_NODES,
                 node_strategy: :best_bound, branch_strategy: :most_fractional)
    types = variable_types || @variable_types
    if types.nil? || types.length != @num_non_slack_vars
      raise ArgumentError,
        "variable_types must list one entry per variable " \
        "(#{@num_non_slack_vars} variables)"
    end
    types.each do |type|
      unless MILP::VARIABLE_TYPES.include?(type)
        raise ArgumentError,
          "unknown variable type #{type.inspect}; " \
          "use :continuous, :integer or :binary"
      end
    end
    unless MILP::NODE_STRATEGIES.include?(node_strategy)
      raise ArgumentError,
        "unknown node_strategy #{node_strategy.inspect}; " \
        "use #{MILP::NODE_STRATEGIES.join(' or ')}"
    end
    unless MILP::BRANCH_STRATEGIES.include?(branch_strategy)
      raise ArgumentError,
        "unknown branch_strategy #{branch_strategy.inspect}; " \
        "use #{MILP::BRANCH_STRATEGIES.join(' or ')}"
    end
    unless max_nodes.is_a?(Integer) && max_nodes >= 1
      raise ArgumentError, "max_nodes must be a positive integer"
    end

    MILP::BranchAndBound.new(
      @objective_coefficients, @constraint_coefficients, @constraint_rhs,
      types,
      max_nodes: max_nodes, node_strategy: node_strategy,
      branch_strategy: branch_strategy
    ).solve
  end
end
