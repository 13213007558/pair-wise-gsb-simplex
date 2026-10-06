require 'matrix'

class Vector
  public :[]=
end

class Simplex
  DEFAULT_MAX_PIVOTS = 10_000

  class UnboundedProblem < StandardError
  end

  class InfeasibleProblem < StandardError
  end

  # A saved optimal basis. Belongs to exactly one model: the fingerprint
  # captures the objective coefficients, the constraint matrix and the
  # problem dimensions, so a basis saved from one model cannot be
  # (mis)applied to another.
  class Basis
    attr_reader :basic_vars, :reduced_costs, :basis_inverse, :fingerprint

    def initialize(basic_vars, reduced_costs, basis_inverse, fingerprint)
      @basic_vars     = basic_vars
      @reduced_costs  = reduced_costs
      @basis_inverse  = basis_inverse
      @fingerprint    = fingerprint
    end
  end

  attr_accessor :max_pivots
  attr_reader :pivot_history

  def initialize(c, a, b)
    @pivot_count = 0
    @max_pivots = DEFAULT_MAX_PIVOTS
    @pivot_history = []

    # Problem dimensions
    @num_non_slack_vars = a.first.length
    @num_constraints    = b.length
    @num_vars           = @num_non_slack_vars + @num_constraints

    # Set up initial matrix A and vectors b, c
    @c = Vector[*c.map {|c1| -1*c1 } + [0]*@num_constraints]
    @a = a.map {|a1| Vector[*(a1.clone + [0]*@num_constraints)]}
    @b = Vector[*b.clone]

    unless @a.all? {|a| a.size == @c.size } and @b.size == @a.length
      raise ArgumentError, "Input arrays have mismatched dimensions" 
    end

    0.upto(@num_constraints - 1) {|i| @a[i][@num_non_slack_vars + i] = 1 }

    # set initial solution: all non-slack variables = 0
    @x          = Vector[*([0]*@num_vars)]
    @basic_vars = (@num_non_slack_vars...@num_vars).to_a
    update_solution

    # Keep pristine copies of the model so a saved basis can be
    # re-expanded into a tableau later (warm start), even after the
    # working tableau has been mutated by pivots.
    @objective_c = c.clone
    @initial_c   = @c.clone
    @initial_a   = @a.map {|row| Vector[*row.to_a] }
    @initial_b   = Vector[*b.clone]
  end

  def solution
    solve
    current_solution
  end

  def current_solution
    @x.to_a[0...@num_non_slack_vars]
  end

  def objective_value
    current_solution.each_with_index.reduce(0) do |sum, (x_i, i)|
      sum + @objective_c[i] * x_i
    end
  end

  def feasible?
    @b.all? {|b_i| b_i >= 0 }
  end

  # Save the current basis (typically the optimal one, after #solve)
  # together with everything needed to rebuild the tableau for a
  # modified right-hand side later.
  def save_basis
    basis_inverse = @a.map {|row| row.to_a[@num_non_slack_vars..-1] }
    Basis.new(@basic_vars.dup, @c.dup, basis_inverse, model_fingerprint)
  end

  # Re-optimize after the right-hand side changed, reusing a previously
  # saved optimal basis. The reduced costs do not depend on b, so the
  # rebuilt tableau is still dual feasible; if some basic variables
  # became negative the dual simplex method restores primal
  # feasibility, then the primal method finishes the job.
  #
  # Raises ArgumentError if the basis was saved from a different model.
  # Raises InfeasibleProblem if the new rhs makes the problem infeasible
  # (the object can still be reused afterwards via another warm_start).
  def warm_start(basis, new_b = nil)
    unless basis.is_a?(Basis) && basis.fingerprint == model_fingerprint
      raise ArgumentError, "Basis was saved from a different model"
    end
    new_b ||= @initial_b.to_a
    unless new_b.length == @num_constraints
      raise ArgumentError, "Right-hand side has wrong dimension"
    end

    @pivot_history = []

    b_inv = Matrix.rows(basis.basis_inverse)
    a0    = Matrix.rows(@initial_a.map {|row| row.to_a })

    @a = (b_inv * a0).row_vectors
    @b = b_inv * Vector[*new_b]
    @c = basis.reduced_costs.dup
    @basic_vars = basis.basic_vars.dup
    update_solution

    dual_solve
    solve
    self
  end

  # Append one more inequality constraint (coeffs . x <= rhs) to the
  # model, then re-optimize. The current optimum stays dual feasible,
  # so if the new constraint is violated the dual simplex method
  # restores feasibility.
  def add_constraint(coeffs, rhs)
    unless coeffs.length == @num_non_slack_vars
      raise ArgumentError, "Constraint has wrong dimension"
    end

    @pivot_history = []

    new_slack_var = @num_vars

    # extend existing rows and the objective with the new slack column
    @c = Vector[*(@c.to_a + [0])]
    @a.each_index {|i| @a[i] = Vector[*(@a[i].to_a + [0])] }

    # express the new row in terms of the current non-basic variables
    new_row = Vector[*(coeffs.clone + [0]*@num_constraints + [1])]
    new_rhs = rhs
    row_indices.each do |row_ix|
      factor = new_row[basic_variable_in_row(row_ix)]
      next if factor == 0
      new_row -= factor * @a[row_ix]
      new_rhs -= factor * @b[row_ix]
    end

    @a << new_row
    @b = Vector[*(@b.to_a + [new_rhs])]
    @x = Vector[*(@x.to_a + [0])]
    @basic_vars << new_slack_var
    @num_constraints += 1
    @num_vars        += 1

    # keep the pristine model in sync so future bases/warm starts
    # refer to the extended model (old bases are invalidated)
    @initial_c = Vector[*(@initial_c.to_a + [0])]
    @initial_a.each_index {|i| @initial_a[i] = Vector[*(@initial_a[i].to_a + [0])] }
    @initial_a << Vector[*(coeffs.clone + [0]*(@num_constraints - 1) + [1])]
    @initial_b = Vector[*(@initial_b.to_a + [rhs])]

    update_solution

    dual_solve
    solve
    self
  end

  def update_solution
    0.upto(@num_vars - 1) {|i| @x[i] = 0 }

    @basic_vars.each do |basic_var|
      row_with_1 = row_indices.detect do |row_ix|
        @a[row_ix][basic_var] == 1
      end
      @x[basic_var] = @b[row_with_1]
    end
  end

  def solve
    while can_improve?
      @pivot_count += 1
      raise "Too many pivots" if @pivot_count > max_pivots 
      pivot
    end
  end

  def can_improve?
    !!entering_variable
  end

  def variables
    (0...@c.size).to_a
  end

  def entering_variable
    variables.select { |var| @c[var] < 0 }.
              min_by { |var| @c[var] }
  end

  def pivot
    pivot_column = entering_variable
    pivot_row    = pivot_row(pivot_column)
    raise UnboundedProblem unless pivot_row
    pivot_at(pivot_row, pivot_column, :primal)
  end

  # Dual simplex: while the current basic solution is primal infeasible
  # (some b < 0) but dual feasible (all reduced costs >= 0), pivot the
  # most violated row out of the basis. Raises InfeasibleProblem when a
  # violated row has no negative coefficient, i.e. the problem itself
  # is infeasible.
  def dual_solve
    while (pivot_row = dual_leaving_row)
      @pivot_count += 1
      raise "Too many pivots" if @pivot_count > max_pivots
      pivot_column = dual_entering_column(pivot_row)
      raise InfeasibleProblem unless pivot_column
      pivot_at(pivot_row, pivot_column, :dual)
    end
  end

  def dual_leaving_row
    row_indices.select {|row_ix| @b[row_ix] < 0 }.
                min_by {|row_ix| @b[row_ix] }
  end

  def dual_entering_column(pivot_row)
    candidates = column_indices.select {|col_ix| @a[pivot_row][col_ix] < 0 }.
                                map {|col_ix| [col_ix, @c[col_ix], @a[pivot_row][col_ix]] }
    column_ix, _, _ = *last_min_by(candidates) {|_, c, a| c / -a }
    column_ix
  end

  def pivot_at(pivot_row, pivot_column, phase)
    leaving_var  = basic_variable_in_row(pivot_row)
    replace_basic_variable(leaving_var => pivot_column)

    pivot_ratio = Rational(1, @a[pivot_row][pivot_column])

    # update pivot row
    @a[pivot_row] *= pivot_ratio
    @b[pivot_row] = pivot_ratio * @b[pivot_row]

    # update objective
    @c -= @c[pivot_column] * @a[pivot_row]

    # update A and B
    (row_indices - [pivot_row]).each do |row_ix|
      r = @a[row_ix][pivot_column]
      @a[row_ix] -= r * @a[pivot_row]
      @b[row_ix] -= r * @b[pivot_row]
    end

    update_solution

    @pivot_history << {
      :phase    => phase,
      :entering => pivot_column,
      :leaving  => leaving_var,
      :row      => pivot_row,
    }
  end

  def replace_basic_variable(hash)
    from, to = hash.keys.first, hash.values.first
    @basic_vars.delete(from)
    @basic_vars << to
    @basic_vars.sort!
  end

  def pivot_row(column_ix)
    row_ix_a_and_b = row_indices.map { |row_ix|
      [row_ix, @a[row_ix][column_ix], @b[row_ix]]
    }.reject { |_, a, b|
      a == 0
    }.reject { |_, a, b|
      (b < 0) ^ (a < 0) # negative sign check
    }
    row_ix, _, _ = *last_min_by(row_ix_a_and_b) { |_, a, b|
      Rational(b, a)
    }
    row_ix
  end

  def basic_variable_in_row(pivot_row)
    column_indices.detect do |column_ix|
      @a[pivot_row][column_ix] == 1 and @basic_vars.include?(column_ix)
    end
  end

  def row_indices
    (0...@a.length).to_a
  end

  def column_indices
    (0...@a.first.size).to_a
  end

  def formatted_tableau
    if can_improve?
      pivot_column = entering_variable
      pivot_row    = pivot_row(pivot_column)
    else
      pivot_row = nil
    end
    num_cols = @c.size + 1
    c = formatted_values(@c.to_a)
    b = formatted_values(@b.to_a)
    a = @a.to_a.map {|ar| formatted_values(ar.to_a) }
    if pivot_row
      a[pivot_row][pivot_column] = "*" + a[pivot_row][pivot_column]
    end
    max = (c + b + a + ["1234567"]).flatten.map(&:size).max
    result = []
    result << c.map {|c| c.rjust(max, " ") }
    a.zip(b) do |arow, brow|
      result << (arow + [brow]).map {|a| a.rjust(max, " ") }
      result.last.insert(arow.length, "|")
    end
    lines = result.map {|b| b.join("  ") }
    max_line_length = lines.map(&:length).max
    lines.insert(1, "-"*max_line_length)
    lines.join("\n")
  end

  def formatted_values(array)
    array.map {|c| "%2.3f" % c }
  end

  # like Enumerable#min_by except if multiple values are minimum 
  # it returns the last
  def last_min_by(array)
    best_element, best_value = nil, nil
    array.each do |element|
      value = yield element
      if !best_element || value <= best_value
        best_element, best_value = element, value
      end
    end
    best_element
  end

  def assert(boolean)
    raise unless boolean
  end

  private

  def model_fingerprint
    [
      @num_non_slack_vars,
      @num_constraints,
      @objective_c,
      @initial_a.map {|row| row.to_a },
    ]
  end

end
