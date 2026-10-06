
simplex
=======

[![Build Status](https://travis-ci.org/danlucraft/simplex.png)](https://travis-ci.org/danlucraft/simplex)

A naive pure-Ruby implementation of the Simplex algorithm for solving linear programming problems. Solves maximizations in standard form.

### Changes

**1.3**: Adds `save_basis`/`warm_start` for reusing an optimal basis when the
right-hand side changes, `update_rhs` and `add_constraint` for modifying the
model, dual-simplex feasibility restoration, `objective_value` and a
`pivot_log` of the re-optimization. Raises `Simplex::InfeasibleProblem` if
the problem is infeasible.

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

### Re-optimization: warm starts and model modifications

When the objective and constraint coefficients stay the same but the
right-hand side changes repeatedly, the optimal basis from the previous
solve can be reused instead of solving from scratch:

    simplex = Simplex.new(c, a, b)
    simplex.solution
    basis = simplex.save_basis

    warm = Simplex.new(c, a, new_b)
    warm.warm_start(basis)  # ArgumentError if the basis belongs to another model
    warm.solution

The same object can also be modified in place:

    simplex.update_rhs(new_b)          # change the right-hand side
    simplex.add_constraint([1, 1], 2)  # append one constraint
    simplex.solution

If the saved basis is no longer feasible for the new right-hand side (or
the added constraint is violated), feasibility is restored with the dual
simplex method. If the problem itself has become infeasible,
`Simplex::InfeasibleProblem` is raised; after fixing the right-hand side
the same object can be solved again.

Every re-optimization records its pivots in `simplex.pivot_log` (entries
with `phase` `:primal` or `:dual`, `entering`/`leaving` variables and the
pivot position), and `simplex.objective_value` returns the optimal value
of the objective function.
