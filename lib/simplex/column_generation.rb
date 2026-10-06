# Column generation on top of Simplex.
#
# Solves covering-type LPs whose variables ("columns") are too numerous
# to enumerate up front:
#
#   min  sum_j cost_j * x_j
#   s.t. sum_j a_ij * x_j >= rhs_i   for each row i
#        x_j >= 0
#
# A restricted master problem (RMP) starts with only a few explicit
# initialization columns. After each solve, the dual prices are handed
# to a pricer, which either returns new columns with negative reduced
# cost or a proven lower bound on the reduced cost of every column not
# yet in the master. Optimality is declared only when no new columns
# are found AND the bound is >= -tolerance.
class ColumnGeneration
  DEFAULT_TOLERANCE      = 1e-9
  DEFAULT_MAX_ITERATIONS = 100

  # The restricted master is infeasible: the current columns cannot
  # cover every constraint. Provide initialization columns (e.g.
  # artificial columns with a high cost) that cover all constraints.
  class InfeasibleMasterError < StandardError; end

  # A column key was reused. Keys give columns a stable identity so
  # that the same column is never appended to the master twice.
  class DuplicateColumnError < StandardError; end

  # The pricer produced no new columns and no reduced-cost bound, so
  # optimality cannot be proven.
  class PricingError < StandardError; end

  # A single column of the master problem.
  #
  # +key+           stable identity; two columns with the same key are
  #                 considered the same column
  # +cost+          objective coefficient
  # +coefficients+  constraint coefficients, one per row of the master
  Column = Struct.new(:key, :cost, :coefficients) do
    # Reduced cost for the minimization master: cost - A^T * duals.
    def reduced_cost(duals)
      total = cost
      coefficients.each_with_index do |a, i|
        total -= a * duals[i]
      end
      total
    end
  end

  # What a pricer returns when called with the current dual prices.
  #
  # +columns+  new columns (typically with negative reduced cost);
  #            may be empty
  # +bound+    proven lower bound on the reduced cost of EVERY column
  #            not yet in the master, or nil if the pricer cannot prove
  #            one (e.g. it ran out of its own budget). Optimality is
  #            only declared when no new columns are found AND
  #            bound >= -tolerance.
  PricingResult = Struct.new(:columns, :bound)

  # Result of #solve.
  #
  # +status+ is :optimal (proven by the reduced-cost bound) or
  # :iteration_limit (out of budget; the solution is feasible for the
  # columns seen so far but NOT proven optimal).
  Result = Struct.new(:status, :objective, :solution, :iterations) do
    def optimal?
      status == :optimal
    end
  end

  attr_reader :columns, :rhs, :tolerance
  attr_accessor :max_iterations

  def initialize(rhs:, pricer:, initial_columns: [],
                 max_iterations: DEFAULT_MAX_ITERATIONS,
                 tolerance: DEFAULT_TOLERANCE)
    @rhs            = rhs
    @pricer         = pricer
    @max_iterations = max_iterations
    @tolerance      = tolerance
    @columns        = []
    @columns_by_key = {}
    initial_columns.each {|column| add_column(column) }
  end

  def add_column(column)
    unless column.coefficients.length == @rhs.length
      raise ArgumentError, "column has #{column.coefficients.length} coefficients, " \
        "expected #{@rhs.length}"
    end
    if @columns_by_key.key?(column.key)
      raise DuplicateColumnError, "a column with key #{column.key.inspect} already exists"
    end
    @columns_by_key[column.key] = column
    @columns << column
    self
  end

  def column?(key)
    @columns_by_key.key?(key)
  end

  def solve
    iterations = 0
    loop do
      duals, primal = solve_restricted_master
      pricing     = @pricer.call(duals)
      new_columns = select_new_columns(pricing, duals)

      if new_columns.empty?
        if pricing.bound && pricing.bound >= -@tolerance
          return Result.new(:optimal, objective_value(primal),
                            solution_hash(primal), iterations)
        end
        raise PricingError, "pricer returned no new columns and no " \
          "reduced-cost bound; optimality is not proven"
      end

      if iterations >= @max_iterations
        # Out of budget: feasible for the columns seen so far, but
        # optimality is not proven.
        return Result.new(:iteration_limit, objective_value(primal),
                          solution_hash(primal), iterations)
      end

      new_columns.each {|column| add_column(column) }
      iterations += 1
    end
  end

  private

  # Columns returned by the pricer that are not yet in the master.
  # Duplicates are dropped (stable identity); a duplicate that still
  # has negative reduced cost means the pricer is misbehaving.
  def select_new_columns(pricing, duals)
    (pricing.columns || []).select do |column|
      if @columns_by_key.key?(column.key)
        if column.reduced_cost(duals) < -@tolerance
          raise DuplicateColumnError, "pricer returned column " \
            "#{column.key.inspect} which is already in the master"
        end
        false
      else
        true
      end
    end
  end

  # The RMP (min c*x, A*x >= rhs) is solved via its dual
  # (max rhs*pi, A^T*pi <= c, pi >= 0), which always has a feasible
  # origin when costs are non-negative. The primal solution is read
  # back from the dual values of that solve.
  def solve_restricted_master
    if @columns.empty?
      raise InfeasibleMasterError, "no columns: provide initialization " \
        "columns that cover every constraint"
    end
    simplex = Simplex.new(
      @rhs.dup,
      @columns.map {|column| column.coefficients.dup },
      @columns.map {|column| column.cost }
    )
    begin
      duals = simplex.solution
    rescue Simplex::UnboundedProblem
      raise InfeasibleMasterError, "restricted master is infeasible: " \
        "the current columns cannot cover every constraint; " \
        "add initialization columns"
    end
    [duals, simplex.dual_values]
  end

  def objective_value(primal)
    total = 0
    @columns.each_with_index do |column, j|
      total += column.cost * primal[j]
    end
    total
  end

  def solution_hash(primal)
    hash = {}
    @columns.each_with_index do |column, j|
      hash[column.key] = primal[j] if primal[j] != 0
    end
    hash
  end
end

