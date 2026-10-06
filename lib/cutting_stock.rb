require 'column_generation'

class CuttingStock
  attr_reader :stock_length, :lengths, :demands

  def initialize(stock_length, lengths, demands)
    @stock_length = stock_length.to_i
    @lengths = lengths.map(&:to_i)
    @demands = demands.map(&:to_i)

    unless @lengths.length == @demands.length && !@lengths.empty?
      raise ArgumentError, "Lengths and demands must be non-empty and have the same size"
    end
    unless @lengths.all? { |length| length > 0 } && @demands.all? { |demand| demand >= 0 }
      raise ArgumentError, "Lengths must be positive and demands must be non-negative"
    end
    unless @lengths.all? { |length| length <= @stock_length }
      raise ArgumentError, "Every requested length must fit in one stock piece"
    end
  end

  def initial_columns
    lengths.each_index.map do |index|
      counts = Array.new(lengths.length, 0)
      counts[index] = 1
      column(counts)
    end
  end

  def solve_dynamic(options = {})
    pricer = DynamicPricer.new(self, options)
    master_options = options.dup
    master_options[:pricer] = pricer
    master = ColumnGeneration.new(demands, initial_columns, master_options)
    [master.solve, pricer]
  end

  def solve_full_enumeration
    columns = all_patterns.map { |counts| column(counts) }
    master = ColumnGeneration.new(demands, columns) do
      ColumnGeneration::PricingDecision.optimal(0)
    end
    [master.solve, columns]
  end

  def all_patterns
    patterns = []
    counts = Array.new(lengths.length, 0)
    enumerate_patterns(0, stock_length, counts, patterns)
    patterns
  end

  def column(counts)
    ColumnGeneration::Column.new(counts, 1, pattern_id(counts))
  end

  def pattern_id(counts)
    "pattern_#{counts.map(&:to_s).join('_')}"
  end

  def used_patterns(result)
    result.solution.map do |id, value|
      column = result.columns.fetch(id)
      next if value <= 0
      [column.coefficients.map(&:to_i), value]
    end.compact.sort_by { |counts, _| counts }
  end

  private

  def enumerate_patterns(index, remaining_length, counts, patterns)
    if index == lengths.length
      patterns << counts.dup unless counts.all?(&:zero?)
      return
    end

    maximum_count = remaining_length / lengths[index]
    0.upto(maximum_count) do |count|
      counts[index] = count
      enumerate_patterns(
        index + 1,
        remaining_length - count * lengths[index],
        counts,
        patterns
      )
    end
    counts[index] = 0
  end

  class DynamicPricer
    attr_reader :calls

    def initialize(problem, options = {})
      @problem = problem
      @call_budget = options[:call_budget]
      @calls = 0
    end

    def call(dual_prices, _master)
      @calls += 1
      if @call_budget && @calls > @call_budget
        return ColumnGeneration::PricingDecision.budget_exceeded
      end

      counts = best_pattern(dual_prices)
      value = dot(counts, dual_prices)
      reduced_cost = Rational(1) - value

      if reduced_cost < 0
        ColumnGeneration::PricingDecision.add(@problem.column(counts))
      else
        ColumnGeneration::PricingDecision.optimal(reduced_cost)
      end
    end

    private

    def best_pattern(dual_prices)
      best_values = Array.new(@problem.stock_length + 1, 0)
      best_counts = Array.new(@problem.stock_length + 1) do
        Array.new(@problem.lengths.length, 0)
      end

      1.upto(@problem.stock_length) do |capacity|
        best_values[capacity] = best_values[capacity - 1]
        best_counts[capacity] = best_counts[capacity - 1].dup

        @problem.lengths.each_with_index do |length, index|
          next if length > capacity

          candidate_counts = best_counts[capacity - length].dup
          candidate_counts[index] += 1
          candidate_value = dot(candidate_counts, dual_prices)
          current_value = best_values[capacity]

          if candidate_value > current_value ||
             (candidate_value == current_value && (candidate_counts <=> best_counts[capacity]) > 0)
            best_values[capacity] = candidate_value
            best_counts[capacity] = candidate_counts
          end
        end
      end

      best_counts.fetch(@problem.stock_length)
    end

    def dot(counts, dual_prices)
      counts.each_with_index.map do |count, index|
        Rational(count) * dual_prices[index].to_r
      end.reduce(0, :+)
    end
  end

  private_constant :DynamicPricer
end
