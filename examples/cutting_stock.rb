#!/usr/bin/env ruby
$LOAD_PATH.unshift(File.expand_path("../../lib", __FILE__))

require 'cutting_stock'

problem = CuttingStock.new(
  100,
  [45, 36, 31],
  [97, 610, 395]
)

dynamic_result, _pricer = problem.solve_dynamic
full_result, all_columns = problem.solve_full_enumeration

unless dynamic_result.objective_value == full_result.objective_value
  raise "Column generation objective does not match full enumeration"
end

puts "Stock length: #{problem.stock_length}"
puts "Requested lengths: #{problem.lengths.inspect}"
puts "Demands: #{problem.demands.inspect}"
puts "Full-enumeration columns: #{all_columns.length}"
puts "Generated master columns: #{dynamic_result.columns.length}"
puts "Pricing calls: #{dynamic_result.pricing_calls}"
puts "Optimal stock pieces: #{dynamic_result.objective_value}"
puts "Patterns used:"
problem.used_patterns(dynamic_result).each do |counts, amount|
  puts "  #{counts.inspect} x #{amount.to_f}"
end
