
simplex
=======

[![Build Status](https://travis-ci.org/danlucraft/simplex.png)](https://travis-ci.org/danlucraft/simplex)

A naive pure-Ruby implementation of the Simplex algorithm for solving linear programming problems. Solves maximizations in standard form.

### Changes

**1.2**: Raises `Simplex::UnboundedProblem` if the problem is unbounded.

### Why?

I wrote this because I needed to solve some small allocation problems for a web game I'm writing, 
and there didn't seem to be any Ruby 2.0 bindings for the "pro" solver libraries, 
and anyway they are hard to install on Heroku.

 * *Use it for*: small LP problems, with feasible origins, when you have trouble loading native or Java solvers,
     and when you can accept not that great performance. 
 * *Don't use it for*: large LP problems, problems with infeasible origins,when you have access to native solvers, when you need very fast solving time.

### Usage

To solve the linear programming problem:

    max x +  y

      2x +  y <= 4
       x + 2y <= 3

       x, y >= 0

Like this:

    > simplex = Simplex.new(
      [1, 1],       # coefficients of objective function
      [             # matrix of inequality coefficients on the lhs ...
        [ 2,  1],
        [ 1,  2],
      ],
      [4, 3]        # .. and the rhs of the inequalities
    )
    > simplex.solution
    => [(5/3), (2/3)]

You can manually iterate the algorithm, and review the tableau at each step. For instance:

    > simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    > puts simplex.formatted_tableau
     -1.000   -1.000    0.000    0.000            
    ----------------------------------------------
     *2.000    1.000    1.000    0.000  |    4.000
      1.000    2.000    0.000    1.000  |    3.000

    > simplex.can_improve?
    => true
    > simplex.pivot
    => [0, 3]

    > puts simplex.formatted_tableau
      0.000   -0.500    0.500    0.000            
    ----------------------------------------------
      1.000    0.500    0.500    0.000  |    2.000
      0.000   *1.500   -0.500    1.000  |    1.000

The asterisk indicates what will be the pivot row and column for the next pivot.

### Column generation

For minimization master problems with equality constraints and non-negative column costs, `ColumnGeneration` keeps only a restricted set of columns and asks a pricer for more columns.

    min  c'x
    s.t. Ax = b
         x >= 0

The pricer receives the current dual prices and returns either new columns with negative reduced cost, or a proven reduced-cost bound.

    require 'column_generation'

    initial_columns = ColumnGeneration.identity_columns(b)

    solver = ColumnGeneration.new(b, initial_columns) do |dual_prices, master|
      column = ColumnGeneration::Column.new(a_j, c_j, "stable-column-id")
      reduced_cost = master.reduced_cost(column, dual_prices)

      if reduced_cost < 0
        ColumnGeneration::PricingDecision.add(column)
      else
        ColumnGeneration::PricingDecision.optimal(reduced_cost)
      end
    end

    result = solver.solve
    result.status
    result.solution
    result.dual_prices
    result.objective_value

Initial columns must be supplied explicitly. They make the otherwise infeasible origin a feasible starting master. Use `ColumnGeneration.identity_columns(b)` for unit-cover columns, or pass domain-specific columns.

Possible terminal statuses are `:optimal`, `:budget_exceeded`, `:pricing_limit`, `:duplicate_column`, and `:no_improving_column`. A budget-exceeded result never reports optimality.

Columns without an explicit id receive a deterministic id generated from their coefficients and cost. Returning the same column again is reported through `last_rejected_column_id` and is not appended twice.

### Cutting-stock example

    require 'cutting_stock'

    problem = CuttingStock.new(100, [45, 36, 31], [97, 610, 395])
    dynamic_result, _pricer = problem.solve_dynamic
    full_result, all_patterns = problem.solve_full_enumeration

`solve_dynamic` starts with one single-item pattern per demanded size. Its pricer solves the unbounded knapsack subproblem with a pure-Ruby dynamic program. `solve_full_enumeration` recursively enumerates every feasible cutting pattern and solves the resulting complete LP as a reference. Both produce the LP relaxation; an integer cutting-stock solution would require additional branching.

Run it with:

    ruby examples/cutting_stock.rb
