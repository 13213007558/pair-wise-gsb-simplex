require 'test/unit'
$LOAD_PATH.unshift(File.expand_path("../../lib", __FILE__))
require 'cutting_stock'

class ColumnGenerationTest < Test::Unit::TestCase
  def test_pricing_generates_column_and_reaches_optimum
    initial_columns = ColumnGeneration.identity_columns([2, 3])
    master = ColumnGeneration.new([2, 3], initial_columns) do |dual_prices, _|
      new_column = ColumnGeneration::Column.new([1, 1], 1, "combined")
      reduced_cost = master.reduced_cost(new_column, dual_prices)
      if reduced_cost < 0
        ColumnGeneration::PricingDecision.add(new_column)
      else
        ColumnGeneration::PricingDecision.optimal(reduced_cost)
      end
    end

    result = master.solve

    assert_equal :optimal, result.status
    assert_equal 3, result.objective_value
    assert_equal [0, 1], result.dual_prices
    assert_equal 2, result.solution.fetch("combined")
    assert_equal 1, result.solution.fetch("initial_1")
  end

  def test_explicit_initial_columns_make_infeasible_origin_solvable
    feasible_without_unit_columns = [
      ColumnGeneration::Column.new([2, 0], 1, "a"),
      ColumnGeneration::Column.new([0, 3], 1, "b")
    ]

    master = ColumnGeneration.new([2, 3], feasible_without_unit_columns) do |dual_prices, _|
      ColumnGeneration::PricingDecision.optimal(
        [1 - 2 * dual_prices[0], 1 - 3 * dual_prices[1]].min
      )
    end

    result = master.solve

    assert result.optimal?
    assert_equal 2, result.objective_value
    assert_equal 1, result.solution.fetch("a")
    assert_equal 1, result.solution.fetch("b")
  end

  def test_budget_exceeded_does_not_announce_optimal
    problem = CuttingStock.new(10, [6, 4], [5, 5])
    result, pricer = problem.solve_dynamic(call_budget: 0)

    assert_equal :budget_exceeded, result.status
    refute result.optimal?
    assert_equal 1, pricer.calls
    assert_equal problem.initial_columns.length, result.columns.length
  end

  def test_budget_exceeded_with_candidate_still_does_not_announce_optimal
    improving_column = ColumnGeneration::Column.new([1, 1], 1, "combined")
    master = ColumnGeneration.new(
      [2, 3],
      ColumnGeneration.identity_columns([2, 3])
    ) do
      ColumnGeneration::PricingDecision.new([improving_column], nil, true)
    end

    result = master.solve

    assert_equal :budget_exceeded, result.status
    refute result.optimal?
    assert_equal 3, result.columns.length
    assert result.solution.key?("combined")
  end

  def test_repeated_column_has_stable_identity_and_is_not_added_twice
    repeated = ColumnGeneration::Column.new([2, 3], 0, "combined")
    initial = ColumnGeneration.identity_columns([2, 3]) + [repeated]
    master = ColumnGeneration.new(
      [2, 3],
      initial,
      max_pricing_calls: 2
    ) do
      ColumnGeneration::PricingDecision.add(
        ColumnGeneration::Column.new([2, 3], 0, "returned_again")
      )
    end

    result = master.solve

    assert_equal :duplicate_column, result.status
    assert_equal "combined", result.last_rejected_column_id
    assert_equal 3, result.columns.length
  end

  def test_proven_reduced_cost_bound_stops_pricing
    master = ColumnGeneration.new(
      [1, 1],
      ColumnGeneration.identity_columns([1, 1])
    ) do
      ColumnGeneration::PricingDecision.optimal(Rational(1, 2))
    end

    result = master.solve

    assert result.optimal?
    assert_equal Rational(1, 2), result.reduced_cost_bound
    assert_equal 1, result.pricing_calls
  end

  def test_initial_columns_are_required
    assert_raise ArgumentError do
      ColumnGeneration.new([1], []) { ColumnGeneration::PricingDecision.optimal(0) }
    end
  end

  def test_bad_initial_columns_raise_infeasibility
    master = ColumnGeneration.new(
      [2],
      [ColumnGeneration::Column.new([0], 1, "cannot_cover_rhs")]
    ) do
      ColumnGeneration::PricingDecision.optimal(0)
    end

    assert_raise ColumnGeneration::InfeasibleInitialMaster do
      master.solve
    end
  end

  def test_cutting_stock_dynamic_programming_matches_full_enumeration
    problem = CuttingStock.new(10, [6, 4, 3], [7, 5, 3])
    dynamic_result, pricer = problem.solve_dynamic
    full_result, all_columns = problem.solve_full_enumeration

    assert dynamic_result.optimal?
    assert_equal full_result.objective_value, dynamic_result.objective_value
    assert full_result.solution.values.all? { |value| value >= 0 }
    assert dynamic_result.columns.length < all_columns.length
    assert pricer.calls <= all_columns.length
  end

  def test_documented_cutting_stock_example_matches_full_enumeration
    problem = CuttingStock.new(100, [45, 36, 31], [97, 610, 395])
    dynamic_result, _pricer = problem.solve_dynamic
    full_result, all_columns = problem.solve_full_enumeration

    assert dynamic_result.optimal?
    assert_equal Rational(1809, 4), dynamic_result.objective_value
    assert_equal full_result.objective_value, dynamic_result.objective_value
    assert_equal 11, all_columns.length
    assert_equal 6, dynamic_result.columns.length
    assert_equal 4, dynamic_result.pricing_calls
  end
end
