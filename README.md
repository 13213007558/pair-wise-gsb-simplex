
simplex
=======

[![Build Status](https://travis-ci.org/danlucraft/simplex.png)](https://travis-ci.org/danlucraft/simplex)

A naive pure-Ruby implementation of the Simplex algorithm for solving linear programming problems. Solves maximizations in standard form.

### Changes

**1.3**: Adds a pure-Ruby Lemke complementary pivot solver for linear complementarity problems.

**1.2**: Raises `Simplex::UnboundedProblem` if the problem is unbounded.

### Linear complementarity

To solve an LCP

    find z >= 0, w = Mz + q >= 0, with z_i * w_i = 0 for every i

use `LCP`:

    > result = LCP.solve(
      [[2, 1], [1, 2]],
      [-1, -1]
    )
    > result.status
    => :solved
    > result.solution
    => [[(1/3), (1/3)], [(0/1), (0/1)]]

The solver uses an exact rational Lemke complementary pivot path rather than a
floating-point penalty objective. `result.iterations` records the entering and
leaving variables, ratios, blocking candidates, basis, and tableau for every
pivot.

`result.status` is one of:

 * `:solved`: either `q >= 0`, or the artificial variable `z0` left the basis;
 * `:secondary_ray`: the driving variable is unblocked. `result.ray` contains the vertex,
   direction, non-negativity evidence, augmented residuals, and complementary products;
 * `:iteration_limit`: the pivot limit was reached.

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
