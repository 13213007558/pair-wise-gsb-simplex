
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

### Lexicographic objectives

`Simplex::Lexicographic` solves several objectives in explicit priority
order instead of guessing one set of "big" weights to add up. Each
objective is optimized only over the solutions that are optimal for every
higher-priority objective, so a lower priority can never make a higher one
worse -- not even by an epsilon. All arithmetic is exact (`Rational`), and
directions can be mixed per level:

    > lex = Simplex::Lexicographic.new(
      [
        { coefficients: [1, 1], direction: :max }, # total output first
        { coefficients: [1, 0], direction: :min }, # then cost
      ],
      [[2, 1], [1, 2]],   # same constraints as above
      [4, 3]
    )
    > lex.solution
    => [(5/3), (2/3)]
    > lex.objective_values
    => [(7/3), (5/3)]

Every level keeps its achieved (exact) objective value and the evidence
that it is optimal given the higher-priority levels:

    > lex.levels[0].value
    => (7/3)
    > lex.levels[0].evidence
    => { reduced_costs: [...], improving_columns: [], optimal: true }

`improving_columns` lists the variables that could still improve that
level without degrading any higher-priority level; at the final solution
it is empty for every level.
