require 'test/unit'
$:.push(File.expand_path("../../lib", __FILE__))
require 'simplex'
require_relative '../examples/cutting_stock'

class ColumnGenerationTest < Test::Unit::TestCase
  Column = ColumnGeneration::Column

  def test_cutting_stock_matches_full_enumeration
    generated = solve_with_column_generation
    enumerated, num_patterns = solve_with_full_enumeration

    assert generated.optimal?
    assert enumerated.optimal?
    assert generated.iterations > 0
    assert_equal enumerated.objective, generated.objective
    assert_equal Rational(54, 7), generated.objective
    assert num_patterns > generated.iterations
  end

  def test_generated_columns_are_a_subset_of_enumerated
    patterns = enumerate_patterns(LENGTHS, STOCK_LENGTH)
    cg = ColumnGeneration.new(
      rhs:             DEMAND,
      pricer:          KnapsackPricer.new(LENGTHS, STOCK_LENGTH),
      initial_columns: initialization_columns(LENGTHS, STOCK_LENGTH)
    )
    cg.solve
    assert cg.columns.length < patterns.length
    cg.columns.each do |column|
      assert patterns.include?(column.coefficients)
    end
  end

  def test_infeasible_master_without_initialization_columns
    cg = ColumnGeneration.new(rhs: [4], pricer: lambda {|_| nil })
    assert_raise ColumnGeneration::InfeasibleMasterError do
      cg.solve
    end
  end

  def test_infeasible_master_when_columns_do_not_cover_a_constraint
    cg = ColumnGeneration.new(
      rhs:             [4, 3],
      pricer:          lambda {|_| nil },
      initial_columns: [Column.new(:only_first, 1, [1, 0])]
    )
    assert_raise ColumnGeneration::InfeasibleMasterError do
      cg.solve
    end
  end

  def test_artificial_initialization_columns_make_master_feasible
    # Artificial identity columns with a high cost are the classic way
    # to start from an otherwise infeasible master.
    big_m = 1_000
    artificial = [
      Column.new([1, 0], big_m, [1, 0]),
      Column.new([0, 1], big_m, [0, 1]),
    ]
    cg = ColumnGeneration.new(
      rhs:             [4, 3],
      pricer:          lambda {|duals|
        ColumnGeneration::PricingResult.new([Column.new([1, 1], 1, [1, 1])], 0)
      },
      initial_columns: artificial
    )
    result = cg.solve
    assert result.optimal?
    assert_equal 4, result.objective
    assert_equal({ [1, 1] => 4 }, result.solution)
  end

  def test_duplicate_column_key_rejected
    cg = ColumnGeneration.new(rhs: [4], pricer: lambda {|_| nil })
    cg.add_column(Column.new(:a, 1, [1]))
    assert_raise ColumnGeneration::DuplicateColumnError do
      cg.add_column(Column.new(:a, 2, [2]))
    end
  end

  def test_pricer_cannot_append_same_column_forever
    # A pricer that keeps offering the same column must not cause an
    # infinite loop: the duplicate is dropped, and with no bound the
    # run fails loudly instead of cycling.
    column = Column.new(:same, 1, [1])
    calls  = 0
    pricer = lambda do |_duals|
      calls += 1
      ColumnGeneration::PricingResult.new([column], nil)
    end
    cg = ColumnGeneration.new(
      rhs:             [4],
      pricer:          pricer,
      initial_columns: [Column.new(:init, 10, [1])]
    )
    assert_raise ColumnGeneration::PricingError do
      cg.solve
    end
    assert_equal 2, calls
  end

  def test_budget_exhausted_is_not_optimal
    # Endless supply of improving columns, no bound: the budget must
    # stop the loop and the result must not claim optimality.
    counter = 0
    pricer  = lambda do |_duals|
      counter += 1
      column = Column.new("c#{counter}", Rational(1, counter + 1), [1])
      ColumnGeneration::PricingResult.new([column], nil)
    end
    cg = ColumnGeneration.new(
      rhs:             [4],
      pricer:          pricer,
      initial_columns: [Column.new(:c0, 4, [1])],
      max_iterations:  3
    )
    result = cg.solve
    assert_equal :iteration_limit, result.status
    assert !result.optimal?
    assert_equal 3, result.iterations
  end

  def test_no_columns_and_no_bound_raises
    pricer = lambda {|_duals| ColumnGeneration::PricingResult.new([], nil) }
    cg = ColumnGeneration.new(
      rhs:             [4],
      pricer:          pricer,
      initial_columns: [Column.new(:init, 10, [1])]
    )
    assert_raise ColumnGeneration::PricingError do
      cg.solve
    end
  end

  def test_bound_proves_optimality_without_new_columns
    pricer = lambda {|_duals| ColumnGeneration::PricingResult.new([], 0) }
    cg = ColumnGeneration.new(
      rhs:             [4],
      pricer:          pricer,
      initial_columns: [Column.new(:only, 2, [1])]
    )
    result = cg.solve
    assert result.optimal?
    assert_equal 8, result.objective
    assert_equal({ only: 4 }, result.solution)
  end

  def test_column_coefficient_length_checked
    cg = ColumnGeneration.new(rhs: [4, 3], pricer: lambda {|_| nil })
    assert_raise ArgumentError do
      cg.add_column(Column.new(:bad, 1, [1]))
    end
  end
end
