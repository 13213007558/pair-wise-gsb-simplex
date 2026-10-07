require 'test/unit'
$:.push(File.expand_path("../../lib", __FILE__))
require 'simplex'

class MILPTest < Test::Unit::TestCase
  def solve(c, a, b, types, **opts)
    Simplex.new(c, a, b).solve_milp(variable_types: types, **opts)
  end

  # Exhaustive enumeration of all integer/binary assignments in the box
  # implied by the variable types; returns the best objective or nil.
  def brute_force(c, a, b, types, domain)
    domains = types.map {|t| t == :binary ? [0, 1] : domain.to_a }
    best = nil
    domains[0].product(*domains[1..-1]).each do |x|
      feasible = a.each_with_index.all? do |row, i|
        row.zip(x).map {|aij, xj| aij * xj }.reduce(0, :+) <= b[i]
      end
      next unless feasible
      obj = c.zip(x).map {|cj, xj| cj * xj }.reduce(0, :+)
      best = obj if best.nil? || obj > best
    end
    best
  end

  def assert_feasible(a, b, x)
    a.each_with_index do |row, i|
      lhs = row.zip(x).map {|aij, xj| aij * xj }.reduce(0, :+)
      assert lhs <= b[i], "constraint #{i} violated by #{x.inspect}"
    end
    assert x.all? {|v| v >= 0 }, "negative variable in #{x.inspect}"
  end

  def assert_integer_feasible(x, types)
    types.each_with_index do |type, j|
      next if type == :continuous
      assert_equal 1, Rational(x[j]).denominator,
        "variable #{j} not integer in #{x.inspect}"
      assert [0, 1].include?(x[j]) if type == :binary
    end
  end

  # The LP relaxation optimum is (5/3, 2/3) with objective 7/3; rounding it
  # gives (2, 1) which violates 2x + y <= 4. The solver must return the
  # true integer optimum instead of a rounded continuous solution.
  def test_integer_optimum_is_not_rounded_relaxation
    a = [[2, 1], [1, 2]]
    b = [4, 3]
    types = [:integer, :integer]
    result = solve([1, 1], a, b, types)
    assert_equal :optimal, result.status
    assert result.optimal?
    assert_equal [1, 1], result.solution
    assert_equal 2, result.objective
    assert_equal 2, result.bound
    assert_equal 3, result.nodes
    assert_feasible(a, b, result.solution)
    assert_integer_feasible(result.solution, types)
  end

  def test_binary_knapsack
    a = [[3, 2, 2]]
    b = [4]
    types = [:binary, :binary, :binary]
    result = solve([5, 4, 3], a, b, types)
    assert_equal :optimal, result.status
    assert_equal [0, 1, 1], result.solution
    assert_equal 7, result.objective
    assert_equal 7, result.bound
    assert_equal 3, result.nodes
    assert_feasible(a, b, result.solution)
    assert_integer_feasible(result.solution, types)
  end

  # x = 1/2 is the only point of the continuous relaxation: feasible as an
  # LP, but no integer solution exists.
  def test_continuous_feasible_but_integer_infeasible
    result = solve([1], [[2], [-2]], [1, -1], [:integer])
    assert_equal :no_integer_solution, result.status
    assert !result.optimal?
    assert_nil result.solution
    assert_nil result.objective
    assert_nil result.bound
    assert_equal 3, result.nodes
  end

  def test_infeasible_relaxation
    result = solve([1], [[1], [-1]], [1, -2], [:integer])
    assert_equal :infeasible, result.status
    assert !result.optimal?
    assert_nil result.solution
    assert_nil result.objective
    assert_nil result.bound
    assert_equal 1, result.nodes
  end

  def test_unbounded_relaxation
    result = solve([1, 1], [[1, -1]], [1], [:integer, :integer])
    assert_equal :unbounded, result.status
    assert !result.optimal?
    assert_nil result.solution
    assert_equal 1, result.nodes
  end

  def test_node_limit
    c = [1, 1]
    a = [[2, 1], [1, 2]]
    b = [4, 3]
    types = [:integer, :integer]

    result = solve(c, a, b, types, max_nodes: 1)
    assert_equal :node_limit, result.status
    assert !result.optimal?
    assert_equal 1, result.nodes
    assert_nil result.solution
    assert_equal Rational(7, 3), result.bound

    # One more node finds an incumbent; the bound stays a valid upper bound.
    result = solve(c, a, b, types, max_nodes: 2)
    assert_equal :node_limit, result.status
    assert !result.optimal?
    assert_equal 2, result.nodes
    assert_equal [1, 1], result.solution
    assert_equal 2, result.objective
    assert_equal Rational(7, 3), result.bound
    assert result.bound >= result.objective
  end

  def test_mixed_integer_and_continuous
    a = [[2, 1], [1, 2]]
    b = [4, 3]
    types = [:integer, :continuous]
    result = solve([0, 1], a, b, types)
    assert_equal :optimal, result.status
    assert_equal [0, Rational(3, 2)], result.solution
    assert_equal Rational(3, 2), result.objective
    assert_equal 1, result.nodes
    assert_feasible(a, b, result.solution)
    assert_integer_feasible(result.solution, types)
  end

  def test_all_continuous_matches_relaxation
    result = solve([1, 1], [[2, 1], [1, 2]], [4, 3],
                   [:continuous, :continuous])
    assert_equal :optimal, result.status
    assert_equal [Rational(5, 3), Rational(2, 3)], result.solution
    assert_equal Rational(7, 3), result.objective
    assert_equal 1, result.nodes
  end

  def test_variable_types_from_constructor
    result = Simplex.new([5, 4, 3], [[3, 2, 2]], [4],
                         variable_types: [:binary] * 3).solve_milp
    assert_equal :optimal, result.status
    assert_equal [0, 1, 1], result.solution
    assert_equal 7, result.objective
  end

  def test_search_strategies_agree_and_are_deterministic
    c = [1, 1]
    a = [[2, 1], [1, 2]]
    b = [4, 3]
    types = [:integer, :integer]

    [:best_bound, :depth_first].each do |node_strategy|
      [:most_fractional, :first].each do |branch_strategy|
        result = solve(c, a, b, types,
                       node_strategy: node_strategy,
                       branch_strategy: branch_strategy)
        assert_equal :optimal, result.status
        assert_equal 2, result.objective
        assert_equal 2, result.bound
        assert result.nodes >= 1
      end
    end

    first  = solve(c, a, b, types)
    second = solve(c, a, b, types)
    assert_equal first.solution, second.solution
    assert_equal first.nodes, second.nodes
  end

  def test_invalid_options
    simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    assert_raise(ArgumentError) { simplex.solve_milp }
    assert_raise(ArgumentError) do
      simplex.solve_milp(variable_types: [:integer])
    end
    assert_raise(ArgumentError) do
      simplex.solve_milp(variable_types: [:integer, :real])
    end
    assert_raise(ArgumentError) do
      simplex.solve_milp(variable_types: [:integer, :integer],
                         node_strategy: :random)
    end
    assert_raise(ArgumentError) do
      simplex.solve_milp(variable_types: [:integer, :integer],
                         branch_strategy: :random)
    end
    assert_raise(ArgumentError) do
      simplex.solve_milp(variable_types: [:integer, :integer], max_nodes: 0)
    end
  end

  def test_brute_force_fixed_models
    models = [
      [[1, 1], [[2, 1], [1, 2]], [4, 3], [:integer, :integer], 0..4],
      [[5, 4, 3], [[3, 2, 2]], [4], [:binary] * 3, 0..1],
      [[3, 2, 5], [[1, 1, 1], [2, 1, 0], [0, 1, 2]], [6, 8, 7],
       [:integer] * 3, 0..6],
      [[2, 3, 1, 4], [[1, 2, 1, 2], [3, 1, 2, 1]], [7, 9],
       [:integer] * 4, 0..7],
      [[4, 3, 5, 2], [[2, 1, 3, 1]], [5], [:binary] * 4, 0..1],
    ]
    models.each do |c, a, b, types, domain|
      result = Simplex.new(c, a, b).solve_milp(variable_types: types)
      expected = brute_force(c, a, b, types, domain)
      if expected.nil?
        assert_equal :no_integer_solution, result.status
      else
        assert_equal :optimal, result.status
        assert_equal expected, result.objective
        assert_equal expected, result.bound
        assert_feasible(a, b, result.solution)
        assert_integer_feasible(result.solution, types)
      end
    end
  end

  def test_brute_force_random_integer_models
    rng = Random.new(42)
    8.times do
      n = 3
      c = Array.new(n) { rng.rand(1..9) }
      a = Array.new(3) { Array.new(n) { rng.rand(1..5) } }
      b = Array.new(3) { rng.rand(5..20) }
      types = [:integer] * n
      result = Simplex.new(c, a, b).solve_milp(variable_types: types)
      expected = brute_force(c, a, b, types, 0..20)
      assert_equal :optimal, result.status
      assert_equal expected, result.objective
      assert_equal expected, result.bound
      assert_feasible(a, b, result.solution)
      assert_integer_feasible(result.solution, types)
    end
  end

  def test_brute_force_random_binary_models
    rng = Random.new(7)
    8.times do
      n = 4
      c = Array.new(n) { rng.rand(1..9) }
      a = Array.new(2) { Array.new(n) { rng.rand(1..5) } }
      b = Array.new(2) { rng.rand(4..10) }
      types = [:binary] * n
      result = Simplex.new(c, a, b).solve_milp(variable_types: types)
      expected = brute_force(c, a, b, types, 0..1)
      assert_equal :optimal, result.status
      assert_equal expected, result.objective
      assert_equal expected, result.bound
      assert_feasible(a, b, result.solution)
      assert_integer_feasible(result.solution, types)
    end
  end
end
