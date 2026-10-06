require 'simplex'

# Solves linear programs with several objectives ordered by strict
# (lexicographic) priority, instead of collapsing them into one
# weighted sum.
#
# Each objective is optimized exactly, in priority order. While a
# lower-priority objective is being optimized, a variable may only
# enter the basis if its reduced cost in every higher-priority
# objective row is exactly zero, so higher-priority objectives can
# never get worse -- not even by a Rational epsilon. All arithmetic
# is exact (inputs are converted with +to_r+).
#
#   lex = LexicographicSimplex.new(
#     [
#       { coefficients: [1, 1, 0], direction: :max }, # first:  max coverage
#       { coefficients: [3, 1, 0], direction: :min }, # then:   min cost
#       { coefficients: [0, 0, 1], direction: :max }, # last:   max fairness
#     ],
#     [[1, 0, 0], [0, 1, 0], [1, 1, 0], [-1, 0, 1], [0, -1, 1]],
#     [6, 6, 10, 0, 0]
#   )
#   lex.solution        # => [4, 6, 4]
#   lex.layers.map(&:value) # => [10, 18, 4]
class LexicographicSimplex < Simplex
  # The recorded outcome of optimizing one priority layer.
  #
  # * +value+         -- exact optimal value of the layer's objective,
  #                      in the objective's own terms (not sign-flipped)
  # * +solution+      -- basic feasible solution when the layer finished
  # * +reduced_costs+ -- the layer's objective row at optimality; the
  #                      evidence that no eligible entering variable can
  #                      improve this layer without degrading earlier ones
  # * +pivots+        -- number of pivots spent on this layer
  Layer = Struct.new(:index, :direction, :coefficients, :value,
                     :solution, :reduced_costs, :pivots, keyword_init: true)

  DIRECTIONS = [:max, :min].freeze

  attr_reader :layers

  # objectives: non-empty array of
  #   { coefficients: [...], direction: :max | :min }   (:direction optional, defaults to :max)
  # in decreasing priority order. +a+ and +b+ are the lhs and rhs of the
  # problem's <= constraints, exactly as for Simplex.
  def initialize(objectives, a, b)
    if objectives.nil? || objectives.empty?
      raise ArgumentError, "At least one objective is required"
    end

    @specs = objectives.map { |objective| normalize_objective(objective) }
    num_objective_vars = @specs.first[:coefficients].length
    unless @specs.all? { |spec| spec[:coefficients].length == num_objective_vars }
      raise ArgumentError, "Objectives have mismatched dimensions"
    end

    @current_layer = 0
    @layers = []

    super(maximized_coefficients(@specs.first),
          a.map { |row| row.map(&:to_r) },
          b.map(&:to_r))

    # One objective row per layer, kept in sync through every pivot.
    # Rows use the Simplex convention: negated coefficients of the
    # maximized objective, so negative entries signal improvement.
    @objective_rows = @specs.map do |spec|
      Vector[*(maximized_coefficients(spec).map { |ci| -ci } + [0] * @num_constraints)]
    end
    @c = @objective_rows[0]
  end

  # Optimizes every layer in priority order and records a Layer for each.
  def solve
    return if @solved

    @specs.each_with_index do |spec, layer_ix|
      @current_layer = layer_ix
      @c = @objective_rows[layer_ix]
      pivots_before = @pivot_count
      super
      record_layer(spec, layer_ix, @pivot_count - pivots_before)
    end
    @solved = true
  end

  # A variable may only enter the basis if it leaves every
  # higher-priority objective exactly unchanged (reduced cost exactly
  # zero); among those, the most improving one for the current layer.
  def entering_variable
    row = @objective_rows[@current_layer]
    variables.
      select { |var| higher_priority_rows.all? { |higher| higher[var] == 0 } }.
      select { |var| row[var] < 0 }.
      min_by { |var| row[var] }
  end

  # Exact optimal value of each objective, in priority order.
  def objective_values
    solve
    @layers.map(&:value)
  end

  def update_objective(pivot_row, pivot_column)
    @objective_rows.each_with_index do |row, ix|
      @objective_rows[ix] = row - row[pivot_column] * @a[pivot_row]
    end
    @c = @objective_rows[@current_layer]
  end

  private

  def normalize_objective(objective)
    direction = objective.fetch(:direction, :max)
    unless DIRECTIONS.include?(direction)
      raise ArgumentError, "Unknown objective direction: #{direction.inspect}"
    end
    { coefficients: objective.fetch(:coefficients).map(&:to_r), direction: direction }
  end

  # The objective this solver actually maximizes for a layer.
  def maximized_coefficients(spec)
    spec[:direction] == :max ? spec[:coefficients] : spec[:coefficients].map { |ci| -ci }
  end

  def higher_priority_rows
    @objective_rows[0...@current_layer]
  end

  def record_layer(spec, layer_ix, pivots)
    solution = current_solution
    value = spec[:coefficients].zip(solution).map { |ci, xi| ci * xi }.inject(0, :+)
    @layers << Layer.new(
      index:         layer_ix,
      direction:     spec[:direction],
      coefficients:  spec[:coefficients],
      value:         value,
      solution:      solution,
      reduced_costs: @objective_rows[layer_ix].to_a,
      pivots:        pivots
    )
  end
end
