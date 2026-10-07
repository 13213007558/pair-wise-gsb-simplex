require 'test/unit'
$:.push(File.expand_path("../../lib", __FILE__))
require 'milp'

class MilpTest < Test::Unit::TestCase
  # Exhaustively enumerates all integer assignments in [0, limit] (or {0,1}
  # for binaries) and returns the best objective value, or nil when no
  # integer feasible point exists. For pure integer/binary models only.
  def brute_force(c, a, b, types, limit = 6)
    ranges = types.map { |t| t == :binary ? [0, 1] : (0..limit).to_a }
    best = nil
    ranges.first.product(*ranges[1..-1]).each do |x|
      feasible = a.each_with_index.all? do |row, i|
        lhs = Rational(0)
        row.each_with_index { |coef, j| lhs += coef.to_r * x[j] }
        lhs <= b[i].to_r
      end
      next unless feasible
      objective = Rational(0)
      c.each_with_index { |coef, j| objective += coef.to_r * x[j] }
      best = objective if best.nil? || objective > best
    end
    best
  end

  def assert_valid_solution(result, c, a, b, types)
    point = result.solution
    assert_not_nil point
    a.each_with_index do |row, i|
      lhs = Rational(0)
      row.each_with_index { |coef, j| lhs += coef.to_r * point[j] }
      assert lhs <= b[i].to_r, "constraint #{i} violated by #{point.inspect}"
    end
    types.each_with_index do |type, j|
      next if type == :continuous
      assert_equal point[j], point[j].floor,
                   "variable #{j} not integral in #{point.inspect}"
      if type == :binary
        assert [0, 1].include?(point[j]),
               "binary variable #{j} out of range in #{point.inspect}"
      end
    end
    objective = Rational(0)
    c.each_with_index { |coef, j| objective += coef.to_r * point[j] }
    assert_equal objective, result.objective.to_r
  end

  # The README example with integer variables. The LP optimum is 7/3 at
  # (5/3, 2/3); rounding it gives (2, 1) which violates 2x + y <= 4, so the
  # solver must actually branch, not round.
  def test_integer_program_beats_rounding
    c = [1, 1]
    a = [[2, 1], [1, 2]]
    b = [4, 3]
    types = [:integer, :integer]

    result = Milp.new(c, a, b, types).solve

    assert_equal :optimal, result.status
    assert result.optimal?
    assert_equal brute_force(c, a, b, types), result.objective
    assert_equal 2, result.objective
    assert_equal result.objective, result.bound
    assert result.nodes > 1, "expected branching, got #{result.nodes} node(s)"
    assert_valid_solution(result, c, a, b, types)
  end

  def test_binary_knapsack
    c = [5, 3, 2]
    a = [[2, 2, 3]]
    b = [4]
    types = [:binary, :binary, :binary]
    result = Milp.new(c, a, b, types).solve

    assert_equal :optimal, result.status
    assert_equal brute_force(c, a, b, types), result.objective
    assert_equal 8, result.objective
    assert_valid_solution(result, c, a, b, types)
  end

  # A binary variable must be capped at 1 even though the constraints
  # would allow x = 2.
  def test_binary_upper_bound_enforced
    result = Milp.new([3], [[1]], [2], [:binary]).solve
    assert_equal :optimal, result.status
    assert_equal [1], result.solution
    assert_equal 3, result.objective
  end

  # Continuously feasible (x = 1/2) but no integer point exists.
  def test_continuously_feasible_but_integer_infeasible
    lp = Milp::LP.solve([1], [[2], [-2]], [1, -1])
    assert_equal :optimal, lp.status
    assert_equal Rational(1, 2), lp.objective

    result = Milp.new([1], [[2], [-2]], [1, -1], [:integer]).solve
    assert_equal :no_integer_solution, result.status
    assert_nil result.solution
    assert_nil result.objective
    assert !result.optimal?
  end

  # The LP relaxation itself is infeasible (x <= 1 and x >= 2).
  def test_lp_infeasible
    result = Milp.new([1], [[1], [-1]], [1, -2], [:integer]).solve
    assert_equal :infeasible, result.status
    assert_nil result.solution
    assert_nil result.bound
  end

  def test_unbounded_relaxation
    result = Milp.new([1, 1], [[-1, 1]], [1], [:integer, :integer]).solve
    assert_equal :unbounded, result.status
  end

  def test_node_limit_exceeded
    c = [1, 1]
    a = [[2, 1], [1, 2]]
    b = [4, 3]
    types = [:integer, :integer]
    result = Milp.new(c, a, b, types, max_nodes: 1).solve

    assert_equal :node_limit_exceeded, result.status
    assert_equal 1, result.nodes
    assert_nil result.solution
    assert !result.optimal?
    # The reported bound is a valid global bound: it dominates the true
    # integer optimum (2) and equals the root LP value 7/3.
    assert_equal Rational(7, 3), result.bound
    assert result.bound >= brute_force(c, a, b, types)
  end

  def test_node_limit_with_incumbent
    c = [3, 2]
    a = [[2, 1], [1, 2]]
    b = [4, 3]
    types = [:integer, :integer]
    full = Milp.new(c, a, b, types).solve
    assert_equal :optimal, full.status

    limited = Milp.new(c, a, b, types, max_nodes: full.nodes - 1).solve
    assert_equal :node_limit_exceeded, limited.status
    assert_equal full.nodes - 1, limited.nodes
    if limited.solution
      assert_valid_solution(limited, c, a, b, types)
      assert limited.objective <= full.objective
    end
    assert limited.bound >= full.objective
  end

  def test_mixed_integer_and_continuous
    # max x + 2y, 2x + 3y <= 7, x integer, y continuous
    result = Milp.new([1, 2], [[2, 3]], [7], [:integer, :continuous]).solve
    assert_equal :optimal, result.status
    assert_equal Rational(14, 3), result.objective
    assert_equal 0, result.solution[0]
    assert_equal Rational(7, 3), result.solution[1]
  end

  def test_all_continuous_reduces_to_lp
    result = Milp.new([1, 1], [[2, 1], [1, 2]], [4, 3]).solve
    assert_equal :optimal, result.status
    assert_equal [Rational(5, 3), Rational(2, 3)], result.solution
    assert_equal Rational(7, 3), result.objective
    assert_equal 1, result.nodes
  end

  # Every combination of the explicit branching/node-selection strategies
  # must solve the same model to the same proven optimum.
  def test_strategies_all_converge
    c = [3, 2, 5]
    a = [[1, 2, 1], [3, 1, 2], [2, 1, 1]]
    b = [9, 12, 8]
    types = [:integer, :integer, :integer]
    expected = brute_force(c, a, b, types)

    [:first, :most_fractional].each do |branch_strategy|
      [:depth_first, :best_bound].each do |node_selection|
        result = Milp.new(c, a, b, types,
                          branch_strategy: branch_strategy,
                          node_selection: node_selection).solve
        assert_equal :optimal, result.status
        assert_equal expected, result.objective,
                     "#{branch_strategy}/#{node_selection}"
        assert_equal result.objective, result.bound
        assert_valid_solution(result, c, a, b, types)
      end
    end
  end

  # Deterministic pseudo-random small models checked against exhaustive
  # enumeration, under both strategy combinations.
  def test_exhaustive_verification_on_small_models
    rng = Random.new(42)
    5.times do |model|
      num_vars = 3
      c = Array.new(num_vars) { rng.rand(1..6) }
      a = Array.new(3) { Array.new(num_vars) { rng.rand(1..4) } }
      b = Array.new(3) { rng.rand(5..12) }
      # explicit upper bounds keep the enumeration finite
      num_vars.times do |j|
        row = Array.new(num_vars) { 0 }
        row[j] = 1
        a << row
        b << 3
      end
      types = [:integer] * num_vars
      expected = brute_force(c, a, b, types, 3)
      assert_not_nil expected, "model #{model} unexpectedly infeasible"

      [:first, :most_fractional].each do |branch_strategy|
        result = Milp.new(c, a, b, types,
                          branch_strategy: branch_strategy).solve
        assert_equal :optimal, result.status, "model #{model}"
        assert_equal expected, result.objective, "model #{model}"
        assert_valid_solution(result, c, a, b, types)
      end
    end
  end

  def test_binary_exhaustive_verification
    rng = Random.new(7)
    3.times do |model|
      num_vars = 4
      c = Array.new(num_vars) { rng.rand(1..9) }
      a = [Array.new(num_vars) { rng.rand(1..5) }]
      b = [rng.rand(6..10)]
      types = [:binary] * num_vars
      expected = brute_force(c, a, b, types)
      assert_not_nil expected, "model #{model} unexpectedly infeasible"

      result = Milp.new(c, a, b, types).solve
      assert_equal :optimal, result.status
      assert_equal expected, result.objective
      assert_valid_solution(result, c, a, b, types)
    end
  end

  def test_invalid_arguments
    assert_raise ArgumentError do
      Milp.new([1, 1], [[1, 1]], [1], [:integer])
    end
    assert_raise ArgumentError do
      Milp.new([1], [[1]], [1], [:bogus])
    end
    assert_raise ArgumentError do
      Milp.new([1], [[1]], [1], [:integer], branch_strategy: :random)
    end
    assert_raise ArgumentError do
      Milp.new([1], [[1]], [1], [:integer], node_selection: :random)
    end
    assert_raise ArgumentError do
      Milp.new([1], [[1]], [1], [:integer], max_nodes: -1)
    end
  end
end
