require 'digest'
require 'simplex'

class ColumnGeneration
  class InfeasibleInitialMaster < StandardError; end
  class PricingError < StandardError; end

  class Column
    attr_reader :id, :coefficients, :cost

    def initialize(coefficients, cost, id = nil)
      @coefficients = coefficients
      @cost = cost
      @id = id
    end

    def with_id(id)
      self.class.new(coefficients, cost, id)
    end
  end

  class PricingDecision
    attr_reader :columns, :reduced_cost_bound, :budget_exceeded

    def initialize(columns, reduced_cost_bound = nil, budget_exceeded = false)
      @columns = Array(columns)
      @reduced_cost_bound = reduced_cost_bound
      @budget_exceeded = budget_exceeded
    end

    def self.add(columns)
      new(Array(columns))
    end

    def self.optimal(reduced_cost_bound = 0)
      if reduced_cost_bound.nil?
        raise ArgumentError, "An optimality decision needs a reduced-cost bound"
      end
      new([], reduced_cost_bound, false)
    end

    def self.budget_exceeded(reduced_cost_bound = nil)
      new([], reduced_cost_bound, true)
    end

    def budget_exceeded?
      @budget_exceeded
    end
  end

  Result = Struct.new(
    :status,
    :columns,
    :solution,
    :dual_prices,
    :objective_value,
    :reduced_cost_bound,
    :pricing_calls,
    :last_rejected_column_id
  ) do
    def optimal?
      status == :optimal
    end
  end

  DEFAULT_TOLERANCE = Rational(1, 1_000_000_000_000)
  DEFAULT_MAX_PRICING_CALLS = 10_000

  attr_reader :rhs, :columns, :pricing_calls

  def initialize(rhs, initial_columns, options = {}, &pricer)
    @rhs = rhs.map(&:to_r).freeze
    @tolerance = Rational(options.fetch(:tolerance, DEFAULT_TOLERANCE).to_r)
    @max_pricing_calls = options.fetch(:max_pricing_calls, DEFAULT_MAX_PRICING_CALLS)
    @pricer = options[:pricer] || pricer
    @pricing_calls = 0
    @columns = []
    @known_column_ids = {}
    @known_signatures = {}

    unless @pricer
      raise ArgumentError, "Supply a pricer object or block"
    end

    initial_columns.each { |column| register_column(column) }
    if @columns.empty?
      raise ArgumentError, "At least one explicit initialization column is required"
    end
  end

  def self.identity_columns(rhs, cost = 1)
    rhs.each_with_index.map do |_, index|
      coefficients = Array.new(rhs.length, 0)
      coefficients[index] = rhs[index].to_r < 0 ? -1 : 1
      Column.new(coefficients, cost, "initial_#{index}")
    end
  end

  def solve
    loop do
      dual_prices, primal_values = solve_restricted_dual
      objective_value = dot(rhs, dual_prices)

      if @max_pricing_calls && @pricing_calls >= @max_pricing_calls
        return result(:pricing_limit, dual_prices, primal_values, objective_value)
      end

      @pricing_calls += 1
      decision = call_pricer(dual_prices)

      if decision.budget_exceeded? && decision.columns.empty?
        return result(
          :budget_exceeded,
          dual_prices,
          primal_values,
          objective_value,
          decision.reduced_cost_bound
        )
      end

      if decision.columns.empty?
        bound = decision.reduced_cost_bound
        if bound && bound >= -@tolerance
          return result(:optimal, dual_prices, primal_values, objective_value, bound)
        end

        return result(
          :no_improving_column,
          dual_prices,
          primal_values,
          objective_value,
          bound
        )
      end

      rejected_id = nil
      accepted = 0

      decision.columns.each do |column|
        duplicate_id = duplicate_column_id(column)
        if duplicate_id
          return result(
            :duplicate_column,
            dual_prices,
            primal_values,
            objective_value,
            decision.reduced_cost_bound,
            duplicate_id
          )
        end

        reduced_cost = reduced_cost(column, dual_prices)
        next unless reduced_cost < -@tolerance

        register_column(column)
        accepted += 1
      end

      if accepted.zero?
        if decision.budget_exceeded?
          return result(
            :budget_exceeded,
            dual_prices,
            primal_values,
            objective_value,
            decision.reduced_cost_bound,
            rejected_id
          )
        end

        return result(
          rejected_id ? :duplicate_column : :no_improving_column,
          dual_prices,
          primal_values,
          objective_value,
          decision.reduced_cost_bound,
          rejected_id
        )
      end

      if decision.budget_exceeded?
        dual_prices, primal_values = solve_restricted_dual
        objective_value = dot(rhs, dual_prices)
        return result(
          :budget_exceeded,
          dual_prices,
          primal_values,
          objective_value,
          decision.reduced_cost_bound,
          rejected_id
        )
      end
    end
  rescue Simplex::UnboundedProblem
    raise InfeasibleInitialMaster,
      "The explicit initialization columns do not make the restricted master feasible"
  end

  def reduced_cost(column, dual_prices = current_dual_prices)
    column.cost.to_r - dot(column.coefficients.map(&:to_r), dual_prices)
  end

  def current_dual_prices
    solve_restricted_dual.first
  end

  def solution
    solve.solution
  end

  private

  class DuplicateColumn < StandardError
    attr_reader :column_id

    def initialize(column_id)
      @column_id = column_id
      super("Column #{column_id} is already in the restricted master")
    end
  end

  def call_pricer(dual_prices)
    decision = @pricer.respond_to?(:call) ?
      @pricer.call(dual_prices, self) :
      @pricer.price(dual_prices, self)

    unless decision.is_a?(PricingDecision)
      raise PricingError, "Pricer must return a ColumnGeneration::PricingDecision"
    end
    decision
  end

  def register_column(column)
    unless column.is_a?(Column)
      raise ArgumentError, "Expected a ColumnGeneration::Column"
    end

    coefficients = column.coefficients.map(&:to_r)
    if coefficients.length != rhs.length
      raise ArgumentError, "Column has #{coefficients.length} coefficients but RHS has #{rhs.length}"
    end

    cost = column.cost.to_r
    if cost < 0
      raise ArgumentError, "Column costs must be non-negative for this solver"
    end
    id = column.id || stable_id(coefficients, cost)
    signature = column_signature(coefficients, cost)

    if @known_signatures.key?(signature)
      raise DuplicateColumn.new(@known_signatures[signature])
    end

    if @known_column_ids.key?(id)
      raise ArgumentError, "Column id #{id.inspect} is already used with different data"
    end

    @known_column_ids[id] = signature
    @known_signatures[signature] = id
    @columns << Column.new(coefficients, cost, id)
    id
  end

  def stable_id(coefficients, cost)
    payload = (coefficients + [cost]).map(&:to_s).join("|")
    "column_#{Digest::SHA1.hexdigest(payload)}"
  end

  def column_signature(coefficients, cost)
    Digest::SHA1.hexdigest((coefficients + [cost]).map(&:to_s).join("|"))
  end

  def duplicate_column_id(column)
    coefficients = column.coefficients.map(&:to_r)
    cost = column.cost.to_r
    signature = column_signature(coefficients, cost)
    if @known_signatures.key?(signature)
      @known_signatures.fetch(signature)
    elsif column.id && @known_column_ids.key?(column.id)
      column.id
    end
  end

  def solve_restricted_dual
    objective = rhs.flat_map { |value| [value, -value] }
    constraints = []

    columns.each do |column|
      constraint = []
      column.coefficients.each do |coefficient|
        constraint.concat([coefficient, -coefficient])
      end
      constraints << constraint
    end

    simplex = Simplex.new(objective, constraints, columns.map(&:cost))
    simplex.solve

    dual_solution = simplex.current_solution
    dual_prices = rhs.each_index.map do |index|
      dual_solution[index * 2] - dual_solution[index * 2 + 1]
    end

    reduced_costs = simplex.objective_row.last(columns.length)
    primal_values = Hash[
      columns.map.with_index do |column, index|
        [column.id, reduced_costs[index]]
      end
    ]

    [dual_prices, primal_values, simplex]
  end

  def result(status, dual_prices, primal_values, objective_value, bound = nil, rejected_id = nil)
    Result.new(
      status,
      Hash[columns.map { |column| [column.id, column] }],
      primal_values,
      dual_prices,
      objective_value,
      bound,
      pricing_calls,
      rejected_id
    )
  end

  def dot(left, right)
    left.each_with_index.map { |value, index| value.to_r * right[index].to_r }.reduce(0, :+)
  end
end
