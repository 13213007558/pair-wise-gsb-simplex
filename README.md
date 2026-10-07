
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

### Mixed-integer and 0/1 variables

`solve_milp` adds branch and bound on top of LP relaxations, so variables can
be declared `:integer` or `:binary` (0/1) in addition to the default
`:continuous`:

    > simplex = Simplex.new(
      [5, 4, 3],       # objective coefficients
      [[3, 2, 2]],     # constraint coefficients
      [4],             # constraint rhs
      variable_types: [:binary, :binary, :binary]
    )
    > result = simplex.solve_milp
    > result.status        # => :optimal
    > result.solution      # => [0, 1, 1]   (best integer solution found)
    > result.objective     # => 7
    > result.bound         # => 7           (global bound on the objective)
    > result.nodes         # => 3           (LP relaxations solved)
    > result.optimal?      # => true        (optimality proven)

The continuous relaxation optimum is never rounded: an integer solution is
only recorded when a relaxation optimum is exactly integer-feasible (all
arithmetic is rational). `result.status` distinguishes the outcomes:

 * `:optimal` - search tree exhausted; `solution` is proven optimal and
   `bound` equals `objective`.
 * `:infeasible` - the LP relaxation itself is infeasible.
 * `:no_integer_solution` - the relaxation is feasible, but the exhausted
   search found no integer-feasible point.
 * `:unbounded` - the LP relaxation is unbounded.
 * `:node_limit` - the node budget ran out; any incumbent is returned in
   `solution` and `bound` remains a valid global bound on the optimum.

The search is deterministic and controlled by explicit strategies:

    simplex.solve_milp(
      max_nodes: 1000,               # node budget (default 10_000)
      node_strategy: :best_bound,    # or :depth_first
      branch_strategy: :most_fractional  # or :first
    )

With `:best_bound` the open node with the largest relaxation bound is
explored first (ties go to the most recently created node, so the down
branch goes first); `:depth_first` always takes the most recently created
node. `:most_fractional` branches on the integer variable whose value is
closest to a half, `:first` on the lowest-indexed fractional one; ties are
broken by variable index.
