require 'test/unit'
$:.push(File.expand_path("../../lib", __FILE__))
require 'lcp'

class LCPTest < Test::Unit::TestCase
  def assert_lcp_solution(result, z, w)
    assert_equal :solved, result.status
    assert_equal z, result.z
    assert_equal w, result.w
    assert((result.z + result.w).all? { |value| value >= 0 })
    assert_equal [Rational(0)] * result.n, result.residuals
    assert_equal [Rational(0)] * result.n, result.augmented_residuals
    assert_equal [Rational(0)] * result.n, result.complementarity_products
    assert_equal Rational(0), result.artificial_value
  end

  def assert_ray_evidence(result)
    ray = result.ray

    assert_equal :secondary_ray, result.status
    assert_equal :driving_variable_unblocked, result.reason
    assert(ray[:non_negative_vertex])
    assert(ray[:non_negative_direction])
    assert(ray[:augmented_residuals].all?(&:zero?))
    assert(ray[:direction_residuals].all?(&:zero?))
    assert(ray[:pivot_column_coefficients].all? { |value| value <= 0 })
    assert(ray[:complementarity].all? { |entry| entry.values_at(:vertex_product, :linear_coefficient, :direction_product).all?(&:zero?) })

    t = Rational(2)
    z = ray[:vertex][:z].zip(ray[:direction][:z]).map { |vertex, direction| vertex + (t * direction) }
    w = ray[:vertex][:w].zip(ray[:direction][:w]).map { |vertex, direction| vertex + (t * direction) }
    artificial = ray[:vertex][:artificial_value] + (t * ray[:direction][:artificial_value])

    assert((z + w + [artificial]).all? { |value| value >= 0 })
    assert(artificial > 0)
    assert_equal [Rational(0)] * result.n, z.zip(w).map { |zi, wi| zi * wi }

    augmented_residuals = result.n.times.map do |i|
      product = result.m[i].each_with_index.inject(Rational(0)) do |sum, (value, j)|
        sum + (value * z[j])
      end
      w[i] - product - artificial - result.q[i]
    end
    assert_equal [Rational(0)] * result.n, augmented_residuals
  end

  def test_initial_feasible_solution
    result = LCP.new([[2, 1], [0, 3]], [2, 4]).solve

    assert_lcp_solution(result, [Rational(0), Rational(0)], [Rational(2), Rational(4)])
    assert_equal :initial_feasible, result.reason
    assert_equal 1, result.iterations.size
    assert_equal :initial_feasible, result.iterations.first[:kind]
  end

  def test_p_matrix_solution_after_artificial_variable_leaves
    result = LCP.new(
      [
        [2, 1],
        [1, 2]
      ],
      [-1, -1]
    ).solve

    third = Rational(1, 3)
    assert_lcp_solution(result, [third, third], [Rational(0), Rational(0)])
    assert_equal :artificial_left_basis, result.reason
    assert result.artificial_exited

    kinds = result.iterations.map { |iteration| iteration[:kind] }
    assert_equal [:artificial, :complementary, :complementary], kinds
    assert_equal [["z0", "w2"], ["z2", "w1"], ["z1", "z0"]],
                 result.iterations.map { |iteration| [iteration.dig(:entering, :name), iteration.dig(:leaving, :name)] }
    assert(result.iterations.all? { |iteration| !iteration.key?(:after) || iteration[:after][:artificial_value] >= 0 })

    final_pivot = result.iterations.last
    assert final_pivot[:artificial_exited]
    assert_equal "z0", final_pivot[:leaving][:name]
  end

  def test_degenerate_zero_ratio_branch_is_deterministic
    result = LCP.new(
      [
        [2, 1],
        [1, 2]
      ],
      [-1, -1]
    )
    result.max_pivots = 10
    result.solve

    artificial_pivot = result.iterations[0]
    degenerate_pivot = result.iterations[1]

    assert_equal 2, artificial_pivot[:candidates].size
    assert_equal "w2", artificial_pivot[:leaving][:name]
    assert_equal Rational(0), degenerate_pivot[:ratio]
    assert_equal [[Rational(1), Rational(-1)], [Rational(0), Rational(-1, 2)]],
                 degenerate_pivot[:candidates].map { |candidate| candidate[:tie_break] }
    assert_equal "w1", degenerate_pivot[:leaving][:name]
  end

  def test_p_matrix_with_zero_initial_slack
    result = LCP.new(
      [
        [1, -1],
        [1,  1]
      ],
      [-1, 0]
    ).solve

    assert_lcp_solution(result, [Rational(1), Rational(0)], [Rational(0), Rational(1)])
    assert_equal :artificial_left_basis, result.reason
    assert_equal [["z0", "w1"], ["z1", "z0"]],
                 result.iterations.select { |iteration| iteration[:kind] != :secondary_ray }.
                                      map { |iteration| [iteration.dig(:entering, :name), iteration.dig(:leaving, :name)] }
  end

  def test_secondary_ray_has_status_and_exact_evidence
    result = LCP.new([[-1]], [-1]).solve

    assert_ray_evidence(result)
    assert_equal "z1", result.ray[:driving_variable][:name]
    assert_equal [Rational(0)], result.z
    assert_equal [Rational(0)], result.w
    assert_equal Rational(1), result.artificial_value
    assert_equal [:artificial, :secondary_ray], result.iterations.map { |iteration| iteration[:kind] }
  end

  def test_degenerate_secondary_ray_keeps_artificial_fixed
    result = LCP.new(
      [
        [0, -1],
        [2,  3]
      ],
      [-3, -4]
    ).solve

    assert_ray_evidence(result)
    assert_equal Rational(0), result.ray[:direction][:artificial_value]
    assert result.ray[:vertex][:artificial_value] > 0
    assert_equal "w2", result.ray[:driving_variable][:name]
    assert(result.ray[:pivot_column_coefficients].all? { |value| value <= 0 })
  end

  def test_iteration_limit_returns_state_without_penalty_objective
    result = LCP.new([[2, 1], [1, 2]], [-1, -1])
    result.max_pivots = 1
    result.solve

    assert_equal :iteration_limit, result.status
    assert_equal :max_pivots_reached, result.reason
    assert_equal 1, result.iterations.count { |iteration| iteration[:kind] == :artificial }
    assert_nil result.ray
  end

  def test_rejects_non_square_matrix
    assert_raise ArgumentError do
      LCP.new([[1, 0]], [1])
    end
  end

  def test_rejects_mismatched_q_length
    assert_raise ArgumentError do
      LCP.new([[1, 0], [0, 1]], [1])
    end
  end
end
