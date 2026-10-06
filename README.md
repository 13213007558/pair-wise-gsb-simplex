
simplex
=======

[![Build Status](https://travis-ci.org/danlucraft/simplex.png)](https://travis-ci.org/danlucraft/simplex)

A naive pure-Ruby implementation of the Simplex algorithm for solving linear programming problems. Solves maximizations in standard form.

### Changes

**1.3**: Adds `Simplex#dual_values` and a `ColumnGeneration` helper for
solving covering-type LPs by column generation (see below).

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

For LPs with more variables ("columns") than you want to enumerate up
front, `ColumnGeneration` solves a covering-type master problem

    min  sum_j cost_j * x_j
    s.t. sum_j a_ij * x_j >= rhs_i
         x_j >= 0

starting from a small set of explicit initialization columns (needed
when the empty or partial column set cannot cover every constraint),
and asks a pricer to find new columns from the dual prices:

    pricer = lambda do |duals|
      column = ... # find a column with negative reduced cost, or prove none exists
      ColumnGeneration::PricingResult.new([column], bound)
    end

    cg = ColumnGeneration.new(
      rhs:             [4],
      pricer:          pricer,
      initial_columns: [ColumnGeneration::Column.new(:init, 1000, [1])]
    )
    result = cg.solve
    result.status    # => :optimal or :iteration_limit
    result.objective # => optimal objective value (when :optimal)
    result.solution  # => { column key => value }

The pricer returns a `PricingResult` with new columns and/or a proven
lower bound on the reduced cost of every column not yet in the master.
Optimality is declared only when no new columns are found AND the bound
is non-negative; if the iteration budget (`max_iterations`) runs out
first, the status is `:iteration_limit` and optimality is not claimed.
Columns have stable identity through their `key`: the same column is
never appended twice.

See `examples/cutting_stock.rb` for a cutting stock solver whose pricer
is an exact dynamic-programming knapsack, compared against a model that
enumerates every cutting pattern up front.
