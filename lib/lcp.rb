class LCP
  DEFAULT_MAX_PIVOTS = 10_000

  attr_accessor :max_pivots
  attr_reader :m, :q, :n, :status, :reason, :iterations, :pivot_count, :z, :w,
              :artificial_value, :artificial_exited, :ray, :message

  def self.solve(m, q)
    lcp = new(m, q)
    lcp.solve
  end

  def initialize(m, q)
    raise ArgumentError, "M must be a non-empty matrix" unless m.is_a?(Array) && !m.empty?
    raise ArgumentError, "q must be an array" unless q.is_a?(Array)

    @n = m.length
    unless m.all? { |row| row.is_a?(Array) && row.length == @n }
      raise ArgumentError, "M must be a square matrix"
    end
    unless q.length == @n
      raise ArgumentError, "q must have the same length as M"
    end

    @m = m.map { |row| row.map { |value| to_rational(value) } }
    @q = q.map { |value| to_rational(value) }
    @max_pivots = DEFAULT_MAX_PIVOTS

    @column_count = (2 * @n) + 1
    @artificial_index = 2 * @n
    @rows = Array.new(@n) { Array.new(@column_count, Rational(0)) }
    @rhs = @q.dup
    @rhs_perturbations = Array.new(@n) { Array.new(@n, Rational(0)) }
    @n.times { |i| @rhs_perturbations[i][i] = Rational(1) }
    @basis = Array.new(@n) { |i| @n + i }
    @pivot_count = 0

    @n.times do |i|
      @n.times { |j| @rows[i][j] = -@m[i][j] }
      @rows[i][@n + i] = Rational(1)
      @rows[i][@artificial_index] = Rational(-1)
    end
  end

  def solve
    return self if @status

    @iterations = []
    @artificial_exited = false
    @ray = nil
    @pivot_count = 0
    @artificial_value = Rational(0)
    update_solution

    if @rhs.all? { |value| value >= 0 }
      @status = :solved
      @reason = :initial_feasible
      @artificial_value = Rational(0)
      @message = "q is non-negative; (z, w) = (0, q) is complementary."
      @iterations << {
        iteration: 0,
        kind: :initial_feasible,
        status: @status,
        reason: @reason,
        after: current_state
      }
      return self
    end

    initial_candidates = pivot_candidates(@artificial_index, true)
    initial_row = select_initial_row(initial_candidates)
    leaving = @basis[initial_row]
    before = current_state

    @iterations << {
      iteration: 1,
      kind: :artificial,
      status: :in_progress,
      entering: variable_info(@artificial_index),
      leaving: variable_info(leaving),
      driving: variable_info(@artificial_index),
      pivot_row: initial_row,
      pivot_column: @artificial_index,
      ratio: @rhs[initial_row] / @rows[initial_row][@artificial_index],
      candidates: initial_candidates,
      before: before
    }
    pivot!(initial_row, @artificial_index)
    update_solution
    @iterations.last[:status] = :in_progress
    @iterations.last[:after] = current_state

    driving = complement(leaving)

    loop do
      if @pivot_count + 1 > @max_pivots
        @status = :iteration_limit
        @reason = :max_pivots_reached
        @message = "Complementary pivoting reached the pivot limit."
        @iterations << {
          iteration: @iterations.size + 1,
          kind: :iteration_limit,
          status: @status,
          driving: variable_info(driving),
          tableau: tableau_snapshot
        }
        break
      end

      candidates = pivot_candidates(driving, false)

      if candidates.empty?
        mark_secondary_ray(driving)
        break
      end

      pivot_row = select_blocking_row(candidates)
      leaving = @basis[pivot_row]
      before = current_state
      iteration_number = @iterations.size + 1

      @iterations << {
        iteration: iteration_number,
        kind: :complementary,
        status: :in_progress,
        entering: variable_info(driving),
        leaving: variable_info(leaving),
        driving: variable_info(driving),
        pivot_row: pivot_row,
        pivot_column: driving,
        ratio: @rhs[pivot_row] / @rows[pivot_row][driving],
        candidates: candidates,
        before: before
      }

      pivot!(pivot_row, driving)
      update_solution

      if leaving == @artificial_index
        @status = :solved
        @reason = :artificial_left_basis
        @artificial_exited = true
        @message = "The artificial variable z0 left the basis at zero."
        @iterations.last[:status] = @status
        @iterations.last[:artificial_exited] = true
      end

      @iterations.last[:after] = current_state
      break if @status == :solved

      driving = complement(leaving)
    end

    self
  end

  def solution
    solve unless @status
    [@z, @w]
  end

  def solved?
    solve unless @status
    @status == :solved
  end

  def secondary_ray?
    solve unless @status
    @status == :secondary_ray
  end

  def residuals
    solve unless @status
    @n.times.map do |i|
      product = @m[i].each_with_index.inject(Rational(0)) do |sum, (value, j)|
        sum + (value * @z[j])
      end
      @w[i] - product - @q[i]
    end
  end

  def augmented_residuals
    solve unless @status
    augmented_lcp_residuals(@z, @w, @artificial_value)
  end

  def complementarity_products
    solve unless @status
    @z.zip(@w).map { |zi, wi| zi * wi }
  end

  def iterations
    solve unless @status
    @iterations
  end

  private

  def to_rational(value)
    unless value.is_a?(Numeric) && !value.is_a?(Complex)
      raise ArgumentError, "Matrix and q entries must be numeric"
    end
    value.to_r
  end

  def variable_name(index)
    if index < @n
      "z#{index + 1}"
    elsif index < @artificial_index
      "w#{index - @n + 1}"
    else
      "z0"
    end
  end

  def variable_info(index)
    if index < @n
      { index: index, name: variable_name(index), pair: index, side: :z }
    elsif index < @artificial_index
      pair = index - @n
      { index: index, name: variable_name(index), pair: pair, side: :w }
    else
      { index: index, name: "z0", pair: nil, side: :artificial }
    end
  end

  def complement(index)
    if index < @n
      index + @n
    elsif index < @artificial_index
      index - @n
    end
  end

  def nonbasic_variables
    (0...@column_count).select { |index| !@basis.include?(index) }
  end

  def pivot_candidates(column, initial_pivot)
    @rows.each_index.select do |row|
      coefficient = @rows[row][column]
      if initial_pivot
        coefficient == -1 && @rhs[row] < 0
      else
        coefficient > 0
      end
    end.map do |row|
      coefficient = @rows[row][column]
      {
        row: row,
        basic_variable: variable_info(@basis[row]),
        coefficient: coefficient,
        ratio: @rhs[row] / coefficient,
        tie_break: @rhs_perturbations[row].map { |value| value / coefficient }
      }
    end
  end

  def select_initial_row(candidates)
    candidates.max_by do |candidate|
      [
        candidate[:ratio],
        *candidate[:tie_break]
      ]
    end[:row]
  end

  def select_blocking_row(candidates)
    candidates.min_by do |candidate|
      [
        candidate[:ratio],
        *candidate[:tie_break]
      ]
    end[:row]
  end

  def pivot!(row, column)
    pivot_value = @rows[row][column]
    @rows[row] = @rows[row].map { |value| value / pivot_value }
    @rhs[row] = @rhs[row] / pivot_value
    @rhs_perturbations[row] = @rhs_perturbations[row].map { |value| value / pivot_value }

    @rows.each_index do |other_row|
      next if other_row == row

      factor = @rows[other_row][column]
      @rows[other_row] = @rows[other_row].each_with_index.map do |value, column_index|
        value - (factor * @rows[row][column_index])
      end
      @rhs[other_row] -= factor * @rhs[row]
      @rhs_perturbations[other_row] = @rhs_perturbations[other_row].each_with_index.map do |value, column_index|
        value - (factor * @rhs_perturbations[row][column_index])
      end
    end

    @basis[row] = column
    @pivot_count += 1
  end

  def update_solution
    values = Array.new(@column_count, Rational(0))
    @basis.each_with_index { |variable, row| values[variable] = @rhs[row] }

    @z = values[0...@n]
    @w = values[@n...(2 * @n)]
    @artificial_value = values[@artificial_index]
  end

  def tableau_snapshot
    {
      form: "basic_variable = rhs - tableau_row.dot(nonbasic_variables)",
      columns: (0...@column_count).map { |index| variable_info(index) },
      rows: @rows.map(&:dup),
      rhs: @rhs.dup,
      rhs_perturbations: @rhs_perturbations.map(&:dup),
      basis: @basis.dup,
      nonbasic: nonbasic_variables
    }
  end

  def current_state
    {
      tableau: tableau_snapshot(),
      basis: @basis.map { |index| variable_info(index) },
      nonbasic: nonbasic_variables.map { |index| variable_info(index) },
      z: @z.dup,
      w: @w.dup,
      artificial_value: @artificial_value
    }
  end

  def mark_secondary_ray(driving)
    vertex = Array.new(@column_count, Rational(0))
    @basis.each_with_index { |variable, row| vertex[variable] = @rhs[row] }

    direction = Array.new(@column_count, Rational(0))
    direction[driving] = Rational(1)
    @basis.each_with_index do |variable, row|
      direction[variable] = -@rows[row][driving]
    end

    vertex_z = vertex[0...@n]
    vertex_w = vertex[@n...(2 * @n)]
    vertex_artificial = vertex[@artificial_index]
    direction_z = direction[0...@n]
    direction_w = direction[@n...(2 * @n)]
    direction_artificial = direction[@artificial_index]

    augmented_residuals = augmented_lcp_residuals(vertex_z, vertex_w, vertex_artificial)
    direction_residuals = homogeneous_augmented_residuals(direction_z, direction_w, direction_artificial)
    complementarity = @n.times.map do |i|
      {
        pair: i,
        vertex_product: vertex_z[i] * vertex_w[i],
        linear_coefficient: (direction_z[i] * vertex_w[i]) + (vertex_z[i] * direction_w[i]),
        direction_product: direction_z[i] * direction_w[i]
      }
    end

    @ray = {
      parameter: "t >= 0",
      driving_variable: variable_info(driving),
      pivot_column_coefficients: @rows.map { |row| row[driving] },
      vertex: {
        z: vertex_z,
        w: vertex_w,
        artificial_value: vertex_artificial
      },
      direction: {
        z: direction_z,
        w: direction_w,
        artificial_value: direction_artificial
      },
      augmented_residuals: augmented_residuals,
      direction_residuals: direction_residuals,
      complementarity: complementarity,
      non_negative_vertex: [*vertex_z, *vertex_w, vertex_artificial].all? { |value| value >= 0 },
      non_negative_direction: [*direction_z, *direction_w, direction_artificial].all? { |value| value >= 0 }
    }

    @status = :secondary_ray
    @reason = :driving_variable_unblocked
    @message = "The driving variable #{variable_name(driving)} is unblocked, giving a secondary ray."
    @iterations << {
      iteration: @iterations.size + 1,
      kind: :secondary_ray,
      status: @status,
      driving: variable_info(driving),
      candidates: [],
      tableau: tableau_snapshot,
      ray: @ray
    }
  end

  def augmented_lcp_residuals(z, w, artificial_value)
    @n.times.map do |i|
      product = @m[i].each_with_index.inject(Rational(0)) do |sum, (value, j)|
        sum + (value * z[j])
      end
      lhs = w[i] - product - artificial_value
      lhs - @q[i]
    end
  end

  def homogeneous_augmented_residuals(z, w, artificial_value)
    @n.times.map do |i|
      product = @m[i].each_with_index.inject(Rational(0)) do |sum, (value, j)|
        sum + (value * z[j])
      end
      w[i] - product - artificial_value
    end
  end
end
