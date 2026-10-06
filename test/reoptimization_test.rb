require 'test/unit'
$:.push(File.expand_path("../../lib", __FILE__))
require 'simplex'

class ReoptimizationTest < Test::Unit::TestCase
  C = [60, 90, 300]
  A = [
    [1, 1, 1],
    [1, 3, 0],
    [2, 0, 1]
  ]
  B = [600, 600, 900]

  def residuals(a, x)
    a.map {|row| row.zip(x).inject(0) {|sum, (aij, xj)| sum + aij * xj } }
  end

  def assert_feasible(a, b, x)
    assert x.all? {|xi| xi >= 0 }, "negative variable in #{x.inspect}"
    residuals(a, x).each_with_index do |lhs, i|
      assert lhs <= b[i], "constraint #{i} violated: #{lhs} > #{b[i]}"
    end
  end

  def test_warm_start_after_rhs_changes_matches_cold_solve
    simplex = Simplex.new(C, A, B)
    simplex.solution
    basis = simplex.save_basis

    [[700, 500, 800], [400, 900, 600], [800, 600, 700]].each do |new_b|
      warm = Simplex.new(C, A, new_b)
      warm.warm_start(basis)
      cold = Simplex.new(C, A, new_b)

      assert_equal cold.objective_value, warm.objective_value
      assert_equal cold.current_solution, warm.current_solution
      assert_equal residuals(A, cold.current_solution),
                   residuals(A, warm.current_solution)
      assert_feasible(A, new_b, warm.current_solution)
    end
  end

  def test_sequential_rhs_updates_on_same_object
    simplex = Simplex.new(C, A, B)
    simplex.solution

    [[700, 500, 800], [800, 600, 700], [400, 900, 600], B].each do |new_b|
      simplex.update_rhs(new_b)
      cold = Simplex.new(C, A, new_b)
      assert_equal cold.objective_value, simplex.objective_value
      assert_feasible(A, new_b, simplex.current_solution)
    end
  end

  def test_warm_start_with_same_rhs_needs_no_pivots
    simplex = Simplex.new(C, A, B)
    simplex.solution
    basis = simplex.save_basis

    warm = Simplex.new(C, A, B)
    warm.warm_start(basis)
    assert_equal simplex.current_solution, warm.solution
    assert_equal [], warm.pivot_log
  end

  def test_dual_simplex_restores_feasibility_after_rhs_change
    simplex = Simplex.new(C, A, B)
    simplex.solution
    basis = simplex.save_basis

    # makes the previous optimal basis primal infeasible but dual feasible
    warm = Simplex.new(C, A, [800, 600, 700])
    warm.warm_start(basis)
    warm.solution

    assert warm.pivot_log.any? {|entry| entry[:phase] == :dual },
           "expected a dual simplex pivot in #{warm.pivot_log.inspect}"
    cold = Simplex.new(C, A, [800, 600, 700])
    assert_equal cold.objective_value, warm.objective_value
    assert_equal cold.current_solution, warm.current_solution
  end

  def test_reoptimization_pivot_log
    simplex = Simplex.new(C, A, B)
    simplex.solution
    simplex.update_rhs([800, 600, 700])
    simplex.solution

    assert_not_empty simplex.pivot_log
    simplex.pivot_log.each do |entry|
      assert_include [:primal, :dual], entry[:phase]
      assert_kind_of Integer, entry[:entering]
      assert_kind_of Integer, entry[:leaving]
      assert_kind_of Integer, entry[:pivot_row]
      assert_kind_of Integer, entry[:pivot_column]
    end
  end

  def test_basis_from_other_model_rejected
    other = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    other.solution

    simplex = Simplex.new(C, A, B)
    assert_raise ArgumentError do
      simplex.warm_start(other.save_basis)
    end
  end

  def test_basis_rejected_after_coefficient_change
    simplex = Simplex.new(C, A, B)
    simplex.solution
    basis = simplex.save_basis

    different_objective = Simplex.new([60, 90, 301], A, B)
    assert_raise ArgumentError do
      different_objective.warm_start(basis)
    end

    different_matrix = Simplex.new(C, [[1, 1, 1], [1, 3, 0], [2, 0, 2]], B)
    assert_raise ArgumentError do
      different_matrix.warm_start(basis)
    end
  end

  def test_invalid_basis_variables_rejected
    simplex = Simplex.new(C, A, B)
    signature = simplex.model_signature
    assert_raise ArgumentError do
      simplex.warm_start(Simplex::Basis.new([3, 4], signature))
    end
    assert_raise ArgumentError do
      simplex.warm_start(Simplex::Basis.new([3, 3, 4], signature))
    end
    assert_raise ArgumentError do
      simplex.warm_start(Simplex::Basis.new([3, 4, 99], signature))
    end
  end

  def test_singular_basis_rejected
    simplex = Simplex.new([1, 1], [[2, 2], [1, 1]], [4, 2])
    bogus = Simplex::Basis.new([0, 1], simplex.model_signature)
    assert_raise Simplex::BasisError do
      simplex.warm_start(bogus)
      simplex.solution
    end
  end

  def test_add_constraint_violated_by_current_solution
    simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    assert_equal [Rational(5, 3), Rational(2, 3)], simplex.solution

    simplex.add_constraint([1, 1], 2)
    simplex.solution

    cold = Simplex.new([1, 1], [[2, 1], [1, 2], [1, 1]], [4, 3, 2])
    assert_equal cold.objective_value, simplex.objective_value
    assert_equal 2, simplex.objective_value
    assert simplex.pivot_log.any? {|entry| entry[:phase] == :dual },
           "expected a dual simplex pivot in #{simplex.pivot_log.inspect}"
    assert_feasible([[2, 1], [1, 2], [1, 1]], [4, 3, 2],
                    simplex.current_solution)
  end

  def test_add_constraint_satisfied_by_current_solution
    simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    simplex.solution
    simplex.add_constraint([1, 1], 10)
    assert_equal [Rational(5, 3), Rational(2, 3)], simplex.solution
    assert_equal [], simplex.pivot_log
  end

  def test_add_constraint_rejects_wrong_dimension
    simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    assert_raise ArgumentError do
      simplex.add_constraint([1, 1, 1], 2)
    end
  end

  def test_basis_saved_before_add_constraint_no_longer_applies
    simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    simplex.solution
    basis = simplex.save_basis
    simplex.add_constraint([1, 1], 2)
    assert_raise ArgumentError do
      simplex.warm_start(basis)
    end
  end

  def test_infeasible_modification_then_feasible_again
    simplex = Simplex.new([1, 1], [[1, 1], [-1, -1]], [2, 0])
    assert_equal 2, simplex.objective_value

    # x + y <= 2 and x + y >= 3 cannot hold simultaneously
    simplex.update_rhs([2, -3])
    assert_raise Simplex::InfeasibleProblem do
      simplex.solution
    end

    # modifying back to a feasible rhs, the same object still solves
    simplex.update_rhs([2, -1])
    assert_equal [2, 0], simplex.solution
    assert_equal 2, simplex.objective_value
    assert_feasible([[1, 1], [-1, -1]], [2, -1], simplex.current_solution)
  end

  def test_update_rhs_rejects_wrong_dimension
    simplex = Simplex.new(C, A, B)
    assert_raise ArgumentError do
      simplex.update_rhs([1, 2])
    end
  end
end
