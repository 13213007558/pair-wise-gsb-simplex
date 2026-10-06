require 'simplex'

class Simplex
  # Solves a sequence of objectives lexicographically, in explicit priority
  # order: each objective is optimized only over the solutions that are
  # optimal for every higher-priority objective. A lower-priority objective
  # can never trade off against a higher-priority one, not even by an
  # epsilon -- which is exactly what weighted-sum scalarization with fixed
  # "big" weights cannot guarantee.
  #
  # All arithmetic is exact (Rational), and every level keeps its achieved
  # objective value plus optimality evidence (the final reduced-cost row and
  # the -- empty -- set of columns that could still improve that level
  # without degrading a higher-priority one).
  #
  # Objectives may freely mix directions:
  #
  #   lex = Simplex::Lexicographic.new([
  #     { coefficients: [1, 0], direction: :min }, # violations first
  #     { coefficients: [0, 1], direction: :min }, # then cost
  #     { coefficients: [1, 1], direction: :max }, # then fairness
  #   ], a, b)
  #   lex.solution         # => [..]        optimal x
  #   lex.objective_values # => [..]        exact value per level
  #   lex.levels           # => [Level, ..] value + evidence per level
  class Lexicographic < Simplex
    Level = Struct.new(:index, :direction, :coefficients, :value, :evidence)

    DIRECTIONS = [:max, :min].freeze

    def initialize(objectives, a, b)
      if objectives.nil? || objectives.empty?
        raise ArgumentError, "at least one objective is required"
      end

      @objectives = objectives.map { |objective| normalize_objective(objective) }
      widths = @objectives.map { |objective| objective[:coefficients].length }.uniq
      if widths.length > 1
        raise ArgumentError, "objectives have mismatched dimensions"
      end

      exact_a = a.map { |row| row.map { |v| exact(v) } }
      exact_b = b.map { |v| exact(v) }

      super(@objectives.first[:row], exact_a, exact_b)

      # Extra objective rows, stored with the same sign convention as @c
      # (negated maximization coefficients, zero slack columns).
      @extra_rows = @objectives[1..-1].map do |objective|
        Vector[*(objective[:row].map { |v| -v } + [0] * @num_constraints)]
      end
      unless @extra_rows.all? { |row| row.size == @c.size }
        raise ArgumentError, "Input arrays have mismatched dimensions"
      end
    end

    # Bland-style rule: the lowest-index column whose reduced-cost vector
    # (highest-priority objective first) is lexicographically negative.
    def entering_variable
      variables.detect { |var| improving_column?(var, @objectives.length - 1) }
    end

    # A column may enter iff the first objective (in priority order) whose
    # reduced cost it changes at all, it changes for the better. Pivoting on
    # such a column therefore never worsens a higher-priority objective in
    # favour of a lower-priority one.
    def improving_column?(var, last_level)
      objective_rows.first(last_level + 1).each do |row|
        coefficient = row[var]
        next if coefficient == 0
        return coefficient < 0
      end
      false
    end

    def update_objective_rows(pivot_row, pivot_column)
      super
      @extra_rows.each_with_index do |row, i|
        @extra_rows[i] = row - row[pivot_column] * @a[pivot_row]
      end
    end

    def levels
      solve
      @levels ||= build_levels
    end

    def objective_values
      levels.map(&:value)
    end

    def objective_rows
      [@c] + @extra_rows
    end

    def build_levels
      x = current_solution
      rows = objective_rows
      @objectives.each_with_index.map do |objective, i|
        signed = dot(objective[:row], x)
        value  = objective[:direction] == :min ? -signed : signed
        improving = variables.select { |var| improving_column?(var, i) }
        Level.new(i, objective[:direction], objective[:coefficients], value, {
          reduced_costs:      rows[i].to_a,
          improving_columns:  improving,
          optimal:            improving.empty?,
        })
      end
    end

    def normalize_objective(objective)
      coefficients, direction =
        case objective
        when Hash
          [objective[:coefficients], objective[:direction] || :max]
        when Array
          [objective[0], objective[1] || :max]
        else
          raise ArgumentError, "objective must be a Hash or a [coefficients, direction] pair"
        end
      unless coefficients.is_a?(Array) && !coefficients.empty?
        raise ArgumentError, "objective is missing coefficients"
      end
      unless DIRECTIONS.include?(direction)
        raise ArgumentError, "direction must be one of #{DIRECTIONS.inspect}"
      end
      coefficients = coefficients.map { |v| exact(v) }
      row = direction == :min ? coefficients.map { |v| -v } : coefficients
      { coefficients: coefficients, direction: direction, row: row }
    end

    def exact(value)
      case value
      when Rational then value
      when Numeric  then value.to_r
      else raise ArgumentError, "non-numeric coefficient: #{value.inspect}"
      end
    end

    def dot(coefficients, x)
      coefficients.zip(x).map { |coeff, val| coeff * val }.reduce(0, :+)
    end
  end
end
