require 'matrix'

class Vector
  public :[]=
end

class Simplex
  DEFAULT_MAX_PIVOTS = 10_000

  class UnboundedProblem < StandardError
  end

  # Raised when the constraints cannot be satisfied simultaneously.
  class InfeasibleProblem < StandardError
  end

  # Raised when a basis cannot be used for the current model.
  class BasisError < ArgumentError
  end

  # A basis saved from a previous solve. Carries the indices of the basic
  # variables plus a signature of the model it belongs to, so it cannot be
  # applied blindly to a different model.
  class Basis
    attr_reader :basic_vars, :signature

    def initialize(basic_vars, signature)
      @basic_vars = basic_vars.freeze
      @signature  = signature
    end
  end

  attr_accessor :max_pivots
  attr_reader :pivot_log, :basic_vars

  def initialize(c, a, b)
    @pivot_count = 0
    @max_pivots = DEFAULT_MAX_PIVOTS

    # Problem dimensions
    @num_non_slack_vars = a.first.length
    @num_constraints    = b.length
    @num_vars           = @num_non_slack_vars + @num_constraints

    # Set up initial matrix A and vectors b, c
    @c = Vector[*c.map {|c1| -1*c1 } + [0]*@num_constraints]
    @a = a.map {|a1| Vector[*(a1.clone + [0]*@num_constraints)]}
    @b = Vector[*b.clone]

    @objective_coeffs = c.clone

    unless @a.all? {|a1| a1.size == @c.size } and @b.size == @a.length
      raise ArgumentError, "Input arrays have mismatched dimensions" 
    end

    0.upto(@num_constraints - 1) {|i| @a[i][@num_non_slack_vars + i] = 1 }

    # Pristine copies of the model, so the tableau can be rebuilt from any
    # basis when warm starting or re-optimizing after a modification.
    @orig_c = Vector[*@c.to_a]
    @orig_a = @a.map {|row| Vector[*row.to_a] }
    @rhs    = Vector[*@b.to_a]

    # set initial solution: all non-slack variables = 0
    @x          = Vector[*([0]*@num_vars)]
    @basic_vars = (@num_non_slack_vars...@num_vars).to_a
    @pivot_log  = []
    @needs_refactorize = false
    update_solution
  end

  def solution
    solve
    current_solution
  end

  def current_solution
    @x.to_a[0...@num_non_slack_vars]
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
    refactorize if @needs_refactorize
    restore_feasibility
    pivot while can_improve?
  end

  def can_improve?
    !!entering_variable
  end

  # Use the dual simplex method to recover primal feasibility when the
  # current basis is dual feasible (all reduced costs non-negative) but
  # the right-hand side has negative entries. This is the typical
  # situation after a warm start with a changed rhs or an added
  # constraint.
  def restore_feasibility
    while (row_ix = dual_leaving_row)
      unless dual_feasible?
        raise BasisError, "basis is neither primal nor dual feasible"
      end
      entering = dual_entering_variable(row_ix)
      raise InfeasibleProblem unless entering
      pivot!(row_ix, entering, :dual)
    end
  end

  def dual_feasible?
    @c.all? {|coeff| coeff >= 0 }
  end

  def primal_feasible?
    @b.all? {|value| value >= 0 }
  end

  def dual_leaving_row
    row_indices.select {|row_ix| @b[row_ix] < 0 }.
                min_by {|row_ix| @b[row_ix] }
  end

  def dual_entering_variable(row_ix)
    column_indices.select {|col_ix| @a[row_ix][col_ix] < 0 }.
                   min_by {|col_ix| @c[col_ix] / -@a[row_ix][col_ix] }
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
    pivot!(pivot_row, pivot_column, :primal)
  end

  def pivot!(pivot_row, pivot_column, phase)
    @pivot_count += 1
    raise "Too many pivots" if @pivot_count > max_pivots 

    leaving_var = basic_variable_in_row(pivot_row)
    @pivot_log << {
      phase:        phase,
      entering:     pivot_column,
      leaving:      leaving_var,
      pivot_row:    pivot_row,
      pivot_column: pivot_column,
    }
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
  end

  # Rebuild the tableau from the pristine model data so that the columns
  # of the current basic variables form an identity matrix.
  def refactorize
    @c = Vector[*@orig_c.to_a]
    @a = @orig_a.map {|row| Vector[*row.to_a] }
    @b = Vector[*@rhs.to_a]

    assigned_rows = []
    @basic_vars.each do |var|
      row_ix = (row_indices - assigned_rows).detect {|ix| @a[ix][var] != 0 }
      raise BasisError, "basis is singular for this model" unless row_ix

      ratio = Rational(1, @a[row_ix][var])
      @a[row_ix] *= ratio
      @b[row_ix] = ratio * @b[row_ix]

      @c -= @c[var] * @a[row_ix]

      (row_indices - [row_ix]).each do |other_ix|
        r = @a[other_ix][var]
        next if r == 0
        @a[other_ix] -= r * @a[row_ix]
        @b[other_ix] -= r * @b[row_ix]
      end

      assigned_rows << row_ix
    end

    @x = Vector[*([0]*@num_vars)]
    @needs_refactorize = false
    update_solution
  end

  # A signature of the model coefficients (objective and constraint
  # matrix, but not the right-hand side). Used to make sure a saved
  # basis is only applied to the model it was saved from.
  def model_signature
    {
      num_non_slack_vars: @num_non_slack_vars,
      num_constraints:    @num_constraints,
      objective:          @orig_c.to_a,
      constraints:        @orig_a.map {|row| row.to_a },
    }
  end

  # Save the current basis (e.g. the optimal basis after solving) so it
  # can be reused later via warm_start.
  def save_basis
    Basis.new(@basic_vars.dup, model_signature)
  end

  # Install a previously saved basis. The basis is only accepted if it
  # was saved from a model with identical coefficients (the right-hand
  # side may differ). The next solve rebuilds the tableau from this
  # basis and, if necessary, uses the dual simplex method to restore
  # feasibility.
  def warm_start(basis)
    raise ArgumentError, "expected a Simplex::Basis" unless basis.is_a?(Basis)
    unless basis.signature == model_signature
      raise ArgumentError, "basis was saved from a different model"
    end
    vars = basis.basic_vars
    unless vars.length == @num_constraints &&
           vars.uniq.length == vars.length &&
           vars.all? {|var| var.is_a?(Integer) && var >= 0 && var < @num_vars }
      raise ArgumentError, "basis variables are invalid for this model"
    end
    @basic_vars = vars.dup
    @pivot_log.clear
    @needs_refactorize = true
    self
  end

  # Replace the right-hand side and re-optimize from the current basis
  # on the next solve.
  def update_rhs(new_rhs)
    unless new_rhs.length == @num_constraints
      raise ArgumentError, "right-hand side has the wrong number of entries"
    end
    @rhs = Vector[*new_rhs.to_a]
    @pivot_log.clear
    @needs_refactorize = true
    self
  end

  # Append one constraint with its own slack variable. If the current
  # solution violates the new constraint, the next solve restores
  # feasibility with the dual simplex method.
  def add_constraint(coeffs, rhs)
    unless coeffs.length == @num_non_slack_vars
      raise ArgumentError, "constraint has the wrong number of coefficients"
    end
    new_var = @num_vars
    @orig_c = Vector[*(@orig_c.to_a + [0])]
    @orig_a.each_with_index do |row, row_ix|
      @orig_a[row_ix] = Vector[*(row.to_a + [0])]
    end
    @orig_a << Vector[*(coeffs.clone + [0]*@num_constraints + [1])]
    @rhs = Vector[*(@rhs.to_a + [rhs])]
    @num_constraints += 1
    @num_vars         += 1
    @basic_vars << new_var
    @x = Vector[*(@x.to_a + [0])]
    @pivot_log.clear
    @needs_refactorize = true
    self
  end

  # Value of the objective function at the optimal solution.
  def objective_value
    solve
    current_solution.zip(@objective_coeffs).inject(0) do |sum, (value, coeff)|
      sum + value * coeff
    end
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

end
