require 'test/unit'
$:.push(File.expand_path("../../lib", __FILE__))
require 'simplex'

class LexicographicTest < Test::Unit::TestCase
  # Exact feasibility check of a solution against the original <= constraints.
  def assert_feasible(a, b, x)
    a.each_with_index do |row, i|
      lhs = row.zip(x).map { |coeff, val| coeff * val }.reduce(0, :+)
      assert lhs <= b[i], "constraint #{i} violated: #{lhs} <= #{b[i]} fails"
    end
    x.each { |val| assert val >= 0, "negative variable in #{x.inspect}" }
  end

  # Scheduling-style trade-off: minimize violations first, then cost.
  #
  # With v = violations and c = cost, the feasible trade-off between the two
  # extreme vertices is:
  #
  #   vertex A: v = 1, c = 0   (one violation, zero cost)
  #   vertex B: v = 0, c = K   (no violations, cost K)
  #
  # Substituting x1 = 1 - v, x2 = K - c turns "v + c/K >= 1" into the
  # standard-form constraint "x1 + x2/K <= 1" with a feasible origin, and
  # "min v, then min c" into "max x1, then max x2".
  #
  # A weighted sum M*v + c scores A as M and B as K, so for EVERY fixed
  # weight M, scaling the cost up to K > M makes the weighted model pick
  # the violating vertex A. The lexicographic model always picks B.
  def weighted_sum_picks(weight, cost_scale)
    a = [[1, Rational(1, cost_scale)]]
    b = [1]
    Simplex.new([weight, 1], a, b).solution
  end

  def lexicographic_picks(cost_scale)
    a = [[1, Rational(1, cost_scale)]]
    b = [1]
    Simplex::Lexicographic.new([
      { coefficients: [1, 0], direction: :max }, # minimize violations
      { coefficients: [0, 1], direction: :max }, # then minimize cost
    ], a, b)
  end

  def test_any_fixed_big_weight_can_be_made_to_choose_wrong
    [100, 10_000, 1_000_000_000].each do |weight|
      cost_scale = 2 * weight # K > M, so the weighted sum flips to vertex A

      weighted = weighted_sum_picks(weight, cost_scale)
      assert_equal [0, cost_scale], weighted,
        "weighted sum with weight #{weight} should pick the violating vertex"

      lex = lexicographic_picks(cost_scale)
      assert_equal [1, 0], lex.solution,
        "lexicographic order must keep zero violations (x1 = 1 - v = 1)"
      assert_equal [1, 0], lex.objective_values
      assert_feasible([[1, Rational(1, cost_scale)]], [1], lex.solution)
    end
  end

  def test_mixed_directions_with_exact_rationals
    a = [
      [1,              0, 0],
      [Rational(1, 3), 0, 1],
      [0,              1, 1],
      [1,              1, 1],
    ]
    b = [2, 4, 4, 6]
    lex = Simplex::Lexicographic.new([
      { coefficients: [1, 0, 0], direction: :max },
      { coefficients: [0, 1, 0], direction: :min },
      { coefficients: [0, 0, 1], direction: :max },
    ], a, b)

    assert_equal [2, 0, Rational(10, 3)], lex.solution
    assert_equal [2, 0, Rational(10, 3)], lex.objective_values
    assert lex.objective_values.all? { |v| v.is_a?(Rational) },
      "objective values must be exact Rationals"
    assert_feasible(a, b, lex.solution)
  end

  def test_levels_keep_exact_values_and_optimality_evidence
    lex = Simplex::Lexicographic.new([
      { coefficients: [1, 0, 0], direction: :max },
      { coefficients: [0, 1, 0], direction: :min },
      { coefficients: [0, 0, 1], direction: :max },
    ], [
      [1,              0, 0],
      [Rational(1, 3), 0, 1],
      [0,              1, 1],
      [1,              1, 1],
    ], [2, 4, 4, 6])

    levels = lex.levels
    assert_equal 3, levels.length
    assert_equal [:max, :min, :max], levels.map(&:direction)
    assert_equal [2, 0, Rational(10, 3)], levels.map(&:value)
    levels.each_with_index do |level, i|
      assert_equal i, level.index
      assert level.value.is_a?(Rational)
      assert_equal [], level.evidence[:improving_columns],
        "level #{i} must have no column that improves it without degrading a higher level"
      assert level.evidence[:optimal]
      assert level.evidence[:reduced_costs].all? { |r| r.is_a?(Rational) }
    end
  end

  def test_lower_priorities_never_regress_higher_ones
    cost_scale = 2_000_000
    a = [[1, Rational(1, cost_scale)]]
    b = [1]

    # Optimum of the top-priority objective solved completely on its own.
    standalone = Simplex.new([1, 0], a, b).solution
    standalone_value = standalone[0]

    lex = lexicographic_picks(cost_scale)
    # Optimizing the second objective afterwards must not move the first
    # one off its standalone optimum, not even by an epsilon.
    assert_equal standalone_value, lex.levels[0].value
    assert_equal 1, lex.levels[0].value
    assert_equal 0, lex.levels[1].value
  end

  def test_accepts_pair_form_and_defaults_to_max
    lex = Simplex::Lexicographic.new([
      [[1, 1]],          # direction defaults to :max
      [[1, 0], :min],
    ], [[2, 1], [1, 2]], [4, 3])
    assert_equal [Rational(5, 3), Rational(2, 3)], lex.solution
    assert_equal :max, lex.levels[0].direction
    assert_equal :min, lex.levels[1].direction
  end

  def test_argument_errors
    assert_raise(ArgumentError) do
      Simplex::Lexicographic.new([], [[1]], [1])
    end
    assert_raise(ArgumentError) do
      Simplex::Lexicographic.new([{ coefficients: [1], direction: :sideways }], [[1]], [1])
    end
    assert_raise(ArgumentError) do
      Simplex::Lexicographic.new([[1, 0], [1]], [[1, 1]], [1])
    end
    assert_raise(ArgumentError) do
      Simplex::Lexicographic.new([[1, 0, 0]], [[1, 1]], [1])
    end
  end
end
