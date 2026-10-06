require 'test/unit'
$:.push(File.expand_path("../../lib", __FILE__))
require 'simplex'

class WarmStartTest < Test::Unit::TestCase
  C = [60, 90, 300]
  A = [
    [1, 1, 1],
    [1, 3, 0],
    [2, 0, 1]
  ]
  B = [600, 600, 900]

  def setup
    @simplex = Simplex.new(C, A, B)
    @simplex.solution
    @basis = @simplex.save_basis
  end

  def fresh(new_b)
    Simplex.new(C, A, new_b)
  end

  def residuals(simplex, a, b)
    x = simplex.current_solution
    a.each_with_index.map do |row, i|
      row.each_with_index.reduce(0) {|sum, (coef, j)| sum + coef * x[j] } - b[i]
    end
  end

  def assert_same_optimum(warm, cold, a, b)
    assert_equal cold.objective_value, warm.objective_value
    residuals(warm, a, b).each do |residual|
      assert residual <= 0, "constraint violated by #{residual}"
    end
  end

  def test_repeated_rhs_changes_reuse_saved_basis
    [[700, 600, 900], [600, 700, 900], [600, 600, 800], [500, 500, 700]].each do |new_b|
      @simplex.warm_start(@basis, new_b)
      assert_same_optimum @simplex, fresh(new_b).tap(&:solve), A, new_b
    end
  end

  def test_warm_start_without_pivots_when_basis_stays_optimal
    @simplex.warm_start(@basis, [700, 600, 900])
    assert_equal [], @simplex.pivot_history
    assert_equal fresh([700, 600, 900]).solution, @simplex.solution
  end

  def test_warm_start_uses_dual_simplex_when_basis_becomes_infeasible
    new_b = [600, 600, 300]
    @simplex.warm_start(@basis, new_b)

    dual_pivots = @simplex.pivot_history.select {|p| p[:phase] == :dual }
    assert dual_pivots.any?, "expected dual simplex pivots, got #{@simplex.pivot_history.inspect}"
    assert_same_optimum @simplex, fresh(new_b).tap(&:solve), A, new_b
  end

  def test_pivot_history_records_entering_and_leaving_variables
    @simplex.warm_start(@basis, [600, 600, 300])
    @simplex.pivot_history.each do |record|
      assert [:primal, :dual].include?(record[:phase])
      assert_kind_of Integer, record[:entering]
      assert_kind_of Integer, record[:leaving]
    end
  end

  def test_basis_from_other_model_is_rejected
    other = Simplex.new([60, 90, 300], [[2, 1, 1], [1, 3, 0], [2, 0, 1]], B)
    other.solution
    foreign_basis = other.save_basis
    assert_raise ArgumentError do
      @simplex.warm_start(foreign_basis, [700, 600, 900])
    end
  end

  def test_basis_from_other_sized_model_is_rejected
    other = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    other.solution
    assert_raise ArgumentError do
      @simplex.warm_start(other.save_basis, B)
    end
  end

  def test_add_constraint_reoptimizes_with_dual_simplex
    simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    assert_equal [Rational(5, 3), Rational(2, 3)], simplex.solution

    simplex.add_constraint([1, 0], 1)
    assert_equal [1, 1], simplex.solution
    assert simplex.pivot_history.any? {|p| p[:phase] == :dual }

    cold = Simplex.new([1, 1], [[2, 1], [1, 2], [1, 0]], [4, 3, 1])
    assert_equal cold.solution, simplex.solution
    assert_equal cold.objective_value, simplex.objective_value
  end

  def test_add_non_binding_constraint_keeps_solution
    simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    solution = simplex.solution
    simplex.add_constraint([1, 1], 100)
    assert_equal solution, simplex.solution
    assert_equal [], simplex.pivot_history
  end

  def test_old_basis_is_invalid_after_add_constraint
    simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    simplex.solution
    stale_basis = simplex.save_basis
    simplex.add_constraint([1, 0], 1)
    assert_raise ArgumentError do
      simplex.warm_start(stale_basis, [4, 3, 2])
    end
    fresh_basis = simplex.save_basis
    simplex.warm_start(fresh_basis, [4, 3, 2])
    cold = Simplex.new([1, 1], [[2, 1], [1, 2], [1, 0]], [4, 3, 2])
    assert_equal cold.solution, simplex.solution
  end

  def test_infeasible_modification_then_recovery
    # max x + y, x + y <= 4, x + y >= 1 (redundant lower bound)
    simplex = Simplex.new([1, 1], [[1, 1], [-1, -1]], [4, 1])
    simplex.solution
    basis = simplex.save_basis

    # x + y <= 4 and x + y >= 5 cannot both hold
    assert_raise Simplex::InfeasibleProblem do
      simplex.warm_start(basis, [4, -5])
    end

    # the same object must still be usable after the failed re-optimization
    simplex.warm_start(basis, [4, -2])
    assert_equal 4, simplex.objective_value
    residuals(simplex, [[1, 1], [-1, -1]], [4, -2]).each do |residual|
      assert residual <= 0, "constraint violated by #{residual}"
    end

    # and a plain re-solve from the original rhs still works too
    simplex.warm_start(basis)
    assert_equal 4, simplex.objective_value
  end

  def test_objective_value_and_feasible_helpers
    assert_equal 180_000, @simplex.objective_value
    assert @simplex.feasible?
  end
end
