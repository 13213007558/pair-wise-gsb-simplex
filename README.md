
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

### Mixed-integer programming

`Milp` (in `lib/milp.rb`) solves mixed-integer linear programs in the same
format, using branch-and-bound on top of an exact two-phase simplex LP
relaxation. Each variable can be `:continuous` (the default), `:integer` or
`:binary`:

    > result = Milp.new(
      [1, 1],              # objective coefficients
      [[2, 1], [1, 2]],    # constraint matrix
      [4, 3],              # right hand side
      [:integer, :integer]
    ).solve
    > result.status
    => :optimal
    > result.solution
    => [2, 0]
    > result.objective
    => 2
    > result.bound
    => 2
    > result.nodes
    => 3

The result is never produced by rounding the continuous optimum; the
returned point is always integer feasible and `bound` is a valid global
bound on the optimal objective. `status` is one of:

 * `:optimal` - the search tree was exhausted and `solution` is proven optimal
   (`bound == objective`)
 * `:no_integer_solution` - the LP relaxation is feasible, but no integer
   feasible point exists
 * `:infeasible` - the LP relaxation itself is infeasible
 * `:unbounded` - the LP relaxation is unbounded
 * `:node_limit_exceeded` - the node budget (`max_nodes:`, default 10,000) was
   exhausted; `solution`/`objective` hold the best integer point found so far
   (if any) and `bound` still holds a valid global bound

The search is deterministic and controlled by explicit strategies:

 * `branch_strategy:` `:most_fractional` (default) or `:first` (fixed
   variable order)
 * `node_selection:` `:best_bound` (default) or `:depth_first`
