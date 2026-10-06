# Cutting stock example for column generation.
#
# Stock rolls of length STOCK_LENGTH must be cut to satisfy demand for
# smaller order lengths. A "column" is a cutting pattern: how many
# pieces of each order length are cut from one stock roll. The LP
# relaxation of the set-covering formulation is:
#
#   min  sum_p x_p
#   s.t. sum_p a_ip * x_p >= demand_i   for each order length i
#        x_p >= 0
#
# Pricing is an unbounded knapsack solved exactly by dynamic
# programming, so the pricer can always prove a reduced-cost bound.
# The generated model is compared against a model in which every
# pattern was enumerated up front.

$:.push(File.expand_path("../../lib", __FILE__))
require 'simplex'

STOCK_LENGTH = 100
ORDERS       = [[45, 7], [36, 5], [31, 4], [14, 8]] # [length, demand]
LENGTHS      = ORDERS.map(&:first)
DEMAND       = ORDERS.map(&:last)

# Pricer: given the dual price of each order length, find the pattern
# with the most negative reduced cost via exact DP over stock
# positions. Because the DP is exact, 1 - best_value is a proven bound
# on the minimum reduced cost, so optimality is always provable.
class KnapsackPricer
  def initialize(lengths, stock_length)
    @lengths      = lengths
    @stock_length = stock_length
  end

  def call(duals)
    best   = Array.new(@stock_length + 1) { 0 }
    choice = Array.new(@stock_length + 1)
    1.upto(@stock_length) do |capacity|
      @lengths.each_with_index do |length, i|
        next if length > capacity
        candidate = best[capacity - length] + duals[i]
        if candidate > best[capacity]
          best[capacity]   = candidate
          choice[capacity] = [i, capacity - length]
        end
      end
    end

    capacity = (0..@stock_length).max_by {|c| best[c] }
    value    = best[capacity]
    counts   = Array.new(@lengths.length) { 0 }
    while capacity > 0 && choice[capacity]
      i, previous = choice[capacity]
      counts[i] += 1
      capacity   = previous
    end

    column = ColumnGeneration::Column.new(counts, 1, counts)
    ColumnGeneration::PricingResult.new([column], 1 - value)
  end
end

# Every pattern fitting into one stock roll, enumerated recursively
# (piece indices are non-decreasing, so each multiset appears once).
def enumerate_patterns(lengths, stock_length)
  patterns  = []
  enumerate = nil
  enumerate = lambda do |from, remaining, counts|
    patterns << counts.dup if counts.any? {|count| count > 0 }
    from.upto(lengths.length - 1) do |i|
      next if lengths[i] > remaining
      counts[i] += 1
      enumerate.call(i, remaining - lengths[i], counts)
      counts[i] -= 1
    end
  end
  enumerate.call(0, stock_length, Array.new(lengths.length) { 0 })
  patterns
end

# Explicit initialization columns: one wasteful pattern per order
# length (cut only that length). They make the restricted master
# feasible from the start.
def initialization_columns(lengths, stock_length)
  lengths.each_with_index.map do |length, i|
    counts = Array.new(lengths.length) { 0 }
    counts[i] = stock_length / length
    ColumnGeneration::Column.new(counts, 1, counts)
  end
end

def solve_with_column_generation
  ColumnGeneration.new(
    rhs:             DEMAND,
    pricer:          KnapsackPricer.new(LENGTHS, STOCK_LENGTH),
    initial_columns: initialization_columns(LENGTHS, STOCK_LENGTH)
  ).solve
end

def solve_with_full_enumeration
  patterns = enumerate_patterns(LENGTHS, STOCK_LENGTH)
  # Null pricer: no columns left to generate, bound 0 proves optimality.
  null_pricer = lambda {|_duals| ColumnGeneration::PricingResult.new([], 0) }
  result = ColumnGeneration.new(
    rhs:             DEMAND,
    pricer:          null_pricer,
    initial_columns: patterns.map {|counts| ColumnGeneration::Column.new(counts, 1, counts) }
  ).solve
  [result, patterns.length]
end

def format_pattern(counts)
  parts = []
  counts.each_with_index do |count, i|
    parts << "#{count}x#{LENGTHS[i]}" if count > 0
  end
  parts.join(" + ")
end

if __FILE__ == $0
  generated = solve_with_column_generation
  enumerated, num_patterns = solve_with_full_enumeration

  puts "orders: #{ORDERS.map {|l, d| "#{d}x length #{l}" }.join(', ')} (stock #{STOCK_LENGTH})"
  puts
  puts "column generation:  status=#{generated.status} " \
       "objective=#{generated.objective.to_f} (#{generated.objective}) " \
       "iterations=#{generated.iterations}"
  generated.solution.each do |counts, value|
    puts "  #{value.to_s.rjust(8)}  #{format_pattern(counts)}"
  end
  puts
  puts "full enumeration:   status=#{enumerated.status} " \
       "objective=#{enumerated.objective.to_f} (#{enumerated.objective}) " \
       "patterns=#{num_patterns}"
  puts
  if generated.objective == enumerated.objective
    puts "objectives match: column generation reached the proven optimum"
  else
    puts "MISMATCH: #{generated.objective} != #{enumerated.objective}"
    exit 1
  end
end

