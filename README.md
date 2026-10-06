
simplex
=======

[![Build Status](https://travis-ci.org/danlucraft/simplex.png)](https://travis-ci.org/danlucraft/simplex)

A naive pure-Ruby implementation of the Simplex algorithm for solving linear programming problems. Solves maximizations in standard form.

### Changes

**1.3**: Warm starts and model updates. `save_basis` captures the optimal basis,
`warm_start(basis, new_b)` re-optimizes after the right-hand side changes
(using the dual simplex method to restore feasibility when needed), and
`add_constraint(coeffs, rhs)` appends a constraint and re-optimizes.
Bases are fingerprinted, so a basis saved from one model cannot be applied
to another. Raises `Simplex::InfeasibleProblem` if a modification makes the
problem infeasible.

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

### Warm starts

When the objective and the constraint matrix stay the same and only the
right-hand side changes repeatedly, reuse the optimal basis instead of
solving from scratch:

    > simplex = Simplex.new([1, 1], [[2, 1], [1, 2]], [4, 3])
    > simplex.solution
    => [(5/3), (2/3)]
    > basis = simplex.save_basis
    > simplex.warm_start(basis, [5, 3])   # re-optimizes for the new rhs
    > simplex.solution
    => [(7/3), (1/3)]
    > simplex.pivot_history               # pivots of the re-optimization
    => []

If the new right-hand side makes the saved basis infeasible, the dual
simplex method restores feasibility automatically (recorded as
`:phase => :dual` entries in `pivot_history`). If the problem itself
became infeasible, `Simplex::InfeasibleProblem` is raised; the object
remains usable and can be warm-started again with a feasible rhs.

A constraint can also be appended to a solved model:

    > simplex.add_constraint([1, 0], 1)   # adds x <= 1 and re-optimizes
    > simplex.solution
    => [1, 1]
