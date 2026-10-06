require 'test/unit'
$:.push(File.expand_path("../../lib", __FILE__))
require 'lcp'

class LCPTest < Test::Unit::TestCase
  # Verifies every LCP condition directly from the definition:
  # z >= 0, w = Mz + q >= 0 and z_i * w_i == 0 for all i.
  def assert_lcp_solution(m, q, result)
    assert result.solved?, "expected a solution, got #{result.status}"
    z, w = result.z, result.w
    assert_equal q.length, z.length
    assert_equal q.length, w.length
    q.length.times do |i|
      wi = Rational(q[i])
      q.length.times {|j| wi += Rational(m[i][j]) * z[j] }
      assert_equal wi, w[i], "w must equal Mz + q"
      assert z[i] >= 0, "z must be nonnegative"
      assert w[i] >= 0, "w must be nonnegative"
      assert_equal 0, z[i] * w[i], "z_i * w_i must be 0"
    end
  end

  def test_trivial_solution_when_q_is_nonnegative
    result = LCP.new([[2, 1], [1, 2]], [3, 4]).solve
    assert result.solved?
    assert_equal [Rational(0), Rational(0)], result.z
    assert_equal [Rational(3), Rational(4)], result.w
    assert_equal [], result.iterations
  end

  def test_p_matrix_2x2
    m = [[2, 1], [1, 2]] # positive definite, hence a P-matrix
    q = [-5, -6]
    result = LCP.new(m, q).solve
    assert_lcp_solution m, q, result
    assert_equal [Rational(4, 3), Rational(7, 3)], result.z
    assert_equal [Rational(0), Rational(0)], result.w
  end

  def test_p_matrix_triangular
    m = [[1, 2], [0, 1]] # positive principal minors: a P-matrix
    q = [-1, -1]
    result = LCP.new(m, q).solve
    assert_lcp_solution m, q, result
    assert_equal [Rational(0), Rational(1)], result.z
    assert_equal [Rational(1), Rational(0)], result.w
  end

  def test_random_p_matrices_always_solved
    srand(801)
    10.times do
      n = 3
      a = Array.new(n) { Array.new(n) { rand(-3..3) } }
      # M = A^T A + nI is symmetric positive definite, hence a P-matrix,
      # so Lemke's algorithm is guaranteed to terminate at a solution.
      m = Array.new(n) do |i|
        Array.new(n) do |j|
          s = (0...n).reduce(0) {|acc, k| acc + a[k][i] * a[k][j] }
          i == j ? s + n : s
        end
      end
      q = Array.new(n) { rand(-10..10) }
      result = LCP.new(m, q).solve
      assert_lcp_solution m, q, result
    end
  end

  def test_degenerate_tied_initial_pivot
    # Both rhs entries are equally negative: the initial pivot row is a
    # tie, and the follow-up ratio test hits a zero (degenerate) rhs.
    m = [[1, 0], [0, 1]]
    q = [-1, -1]
    result = LCP.new(m, q).solve
    assert_lcp_solution m, q, result
    assert_equal [Rational(1), Rational(1)], result.z
    assert_equal [Rational(0), Rational(0)], result.w
    # the path really went through a degenerate (zero-ratio) pivot
    zero_ratios = result.iterations.flat_map {|it| it[:ratios] || [] }.
                  select {|r| r[:ratio] == 0 }
    assert !zero_ratios.empty?, "expected a degenerate pivot in the path"
  end

  def test_degenerate_zero_rhs
    m = [[1, 0], [0, 1]]
    q = [0, -1]
    result = LCP.new(m, q).solve
    assert_lcp_solution m, q, result
    assert_equal [Rational(0), Rational(1)], result.z
    assert_equal [Rational(0), Rational(0)], result.w
  end

  def test_degenerate_diagonal_tie
    m = [[2, 0], [0, 1]]
    q = [-2, -2]
    result = LCP.new(m, q).solve
    assert_lcp_solution m, q, result
    assert_equal [Rational(1), Rational(2)], result.z
    assert_equal [Rational(0), Rational(0)], result.w
  end

  def test_secondary_ray_for_zero_matrix
    # w = q < 0 can never be feasible; Lemke stops on a secondary ray.
    result = LCP.new([[0, 0], [0, 0]], [-1, -2]).solve
    assert_equal :secondary_ray, result.status
    assert !result.solved?
    assert_nil result.z
    ray = result.ray
    assert_equal "z2", ray[:entering]
    assert ray[:column].all? {|c| c <= 0 }, "ray column must have no positive entry"
    assert_equal 2, ray[:column].length
    assert_equal 2, ray[:basic_variables].length
    assert_equal ray[:basic_variables].length, ray[:basic_values].length
    assert ray[:message].is_a?(String)
    assert !result.iterations.empty?
  end

  def test_secondary_ray_for_negative_diagonal
    result = LCP.new([[-1, 0], [0, -1]], [-1, -1]).solve
    assert_equal :secondary_ray, result.status
    assert result.ray[:column].all? {|c| c <= 0 }
  end

  def test_iteration_records_explain_the_pivot_path
    result = LCP.new([[2, 1], [1, 2]], [-5, -6]).solve
    assert result.solved?
    iterations = result.iterations
    assert_equal 3, iterations.size

    first = iterations.first
    assert_equal 1, first[:step]
    assert_equal "z0", first[:entering]
    assert_equal "w2", first[:leaving] # most negative rhs row
    assert_match(/initial pivot/, first[:note])

    # every later entering variable is the complement of the last leaver
    iterations.each_cons(2) do |prev, nxt|
      assert_equal "#{prev[:leaving][0] == 'w' ? 'z' : 'w'}#{prev[:leaving][1..-1]}",
                   nxt[:entering]
    end

    last = iterations.last
    assert_equal "z0", last[:leaving], "solution requires the artificial variable to leave"
    assert_match(/solution found/, last[:note])

    iterations.each do |it|
      assert it[:basis].is_a?(Array)
      assert_equal 2, it[:basic_values].length
      assert it[:note].is_a?(String)
    end
  end

  def test_class_method_solve
    result = LCP.solve([[1, 2], [0, 1]], [-1, -1])
    assert result.solved?
    assert_equal [Rational(0), Rational(1)], result.z
  end

  def test_error_on_mismatched_dimensions
    assert_raise ArgumentError do
      LCP.new([[1, 2]], [1, 2])
    end
    assert_raise ArgumentError do
      LCP.new([[1, 2], [3, 4]], [1])
    end
    assert_raise ArgumentError do
      LCP.new([[1], [2]], [1, 2])
    end
  end

  def test_formatted_tableau_is_a_string
    lcp = LCP.new([[2, 1], [1, 2]], [-5, -6])
    assert lcp.formatted_tableau.is_a?(String)
    lcp.solve
    assert lcp.formatted_tableau.is_a?(String)
  end
end
