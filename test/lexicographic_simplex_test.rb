require 'test/unit'
$:.push(File.expand_path("../../lib", __FILE__))
require 'simplex'
require 'lexicographic_simplex'

class LexicographicSimplexTest < Test::Unit::TestCase
  def assert_constraints_satisfied(a, b, solution)
    a.zip(b).each_with_index do |(row, rhs), ix|
      lhs = row.zip(solution).map { |ai, xi| ai * xi }.inject(0, :+)
      assert lhs <= rhs, "constraint violated: #{lhs} <= #{rhs}"
    end
    assert solution.all? { |xi| xi >= 0 }, "negative variable in #{solution}"
  end

  def objective_value(coefficients, solution)
    coefficients.zip(solution).map { |ci, xi| ci * xi }.inject(0, :+)
  end

  def test_priority_order_changes_the_answer
    a = [[1, 1], [1, 0], [0, 1]]
    b = [4, 2, 3]
    x_first = LexicographicSimplex.new(
      [{ coefficients: [1, 0] }, { coefficients: [0, 1] }], a, b)
    assert_equal [2, 2], x_first.solution
    assert_equal [2, 2], x_first.objective_values
    y_first = LexicographicSimplex.new(
      [{ coefficients: [0, 1] }, { coefficients: [1, 0] }], a, b)
    assert_equal [1, 3], y_first.solution
    assert_equal [3, 1], y_first.objective_values
  end

  # Scheduling-flavoured model: first maximize coverage (i.e. minimize
  # uncovered demand), then minimize cost, then maximize a fairness
  # level bounded by both crews' loads. Directions are mixed.
  def test_mixed_directions_scheduling
    # variables: a, b (work assigned to two crews), f (fairness level)
    a = [
      [ 1,  0, 0],  # a <= 6
      [ 0,  1, 0],  # b <= 6
      [ 1,  1, 0],  # a + b <= 10 (demand)
      [-1,  0, 1],  # f <= a
      [ 0, -1, 1],  # f <= b
    ]
    b = [6, 6, 10, 0, 0]
    objectives = [
      { coefficients: [1, 1, 0], direction: :max }, # coverage
      { coefficients: [3, 1, 0], direction: :min }, # cost
      { coefficients: [0, 0, 1], direction: :max }, # fairness
    ]

    lex = LexicographicSimplex.new(objectives, a, b)
    assert_equal [4, 6, 4], lex.solution
    assert_equal [10, 18, 4], lex.objective_values

    # every layer did real work
    assert_equal [2, 1, 1], lex.layers.map(&:pivots)

    # each layer kept its exact optimal value: later layers did not
    # degrade earlier ones at all
    final = lex.solution
    lex.layers.each do |layer|
      assert_equal layer.value, objective_value(layer.coefficients, final)
    end

    # optimality evidence is recorded per layer
    lex.layers.each_with_index do |layer, ix|
      assert_equal ix, layer.index
      assert_equal objectives[ix][:direction], layer.direction
      assert_equal 3 + 5, layer.reduced_costs.size
    end

    assert_constraints_satisfied(a, b, final)
  end

  # For ANY fixed "big" weight M there is a scale at which the weighted
  # sum M*x + y trades away the top-priority objective for a gain in the
  # second one. Lexicographic optimization never does, at any scale,
  # and keeps the difference exact (Rational).
  def test_any_fixed_big_weight_can_choose_wrong
    [100, 10_000, 1_000_000].each do |m|
      eps = Rational(1, m)
      # vertices: (0,0), (1,0), (1 - 1/m, m), (0, m)
      a = [
        [1, 0],   # x <= 1
        [0, 1],   # y <= m
        [m, eps], # m*x + y/m <= m
      ]
      b = [1, m, m]
      objectives = [
        { coefficients: [1, 0], direction: :max },
        { coefficients: [0, 1], direction: :max },
      ]

      lex = LexicographicSimplex.new(objectives, a, b)
      assert_equal [1, 0], lex.solution
      assert_equal [1, 0], lex.objective_values
      assert_constraints_satisfied(a, b, lex.solution)

      # the weighted sum sacrifices exactly 1/m of the top objective
      weighted = Simplex.new([m, 1], a, b).solution
      assert_equal [Rational(m - 1, m), m], weighted
      assert_equal Rational(1, m), lex.solution[0] - weighted[0]
    end
  end

  def test_rational_coefficients_stay_exact
    lex = LexicographicSimplex.new(
      [
        { coefficients: [Rational(1, 3), Rational(2, 3)], direction: :max },
        { coefficients: [1, -1], direction: :min },
      ],
      [[1, 2], [4, 1]],
      [4, 7]
    )
    solution = lex.solution
    assert solution.all? { |xi| xi.is_a?(Rational) || xi.is_a?(Integer) }
    # layer 1 optimum: x/3 + 2y/3 = 4/3 along the whole edge from
    # (0, 2) to (10/7, 9/7); layer 2 then picks (0, 2) with x - y = -2
    assert_equal [Rational(4, 3), -2], lex.objective_values
    assert_equal [0, 2], solution
  end

  def test_unbounded_layer_raises
    assert_raise Simplex::UnboundedProblem do
      LexicographicSimplex.new(
        [{ coefficients: [1, 0] }, { coefficients: [0, 1] }],
        [[1, 0]],
        [1]
      ).solution
    end
  end

  def test_errors
    assert_raise ArgumentError do
      LexicographicSimplex.new([], [[1]], [1])
    end
    assert_raise ArgumentError do
      LexicographicSimplex.new(
        [{ coefficients: [1, 0], direction: :sideways }], [[1, 0]], [1])
    end
    assert_raise ArgumentError do
      LexicographicSimplex.new(
        [{ coefficients: [1, 0] }, { coefficients: [1, 0, 0] }], [[1, 0]], [1])
    end
  end
end
