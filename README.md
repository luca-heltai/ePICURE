ePICURE
======================================================================

Python Isogeometric CUrve REconstruction
----------------------------------------------------------------------

This is the final project for the course **Applied Mathematics: an
Introduction to Numerical Analysis and Sientific Computing**. It
consists of several small sub-projects, that will combine together in
an open source software for curve resconstruction, analysis and
simulation. We will apply this to micro-swimming problems, but several
other applications are possible...

The students should be familiar with

- Interpolation
- Integration
- Direct and iterative solution of linear systems
- Solution of non-linear systems
- Systems of Ordinary Differential Equations (ODEs)
- Finite difference and Galerkin methods
- Git

Objects
------------------------------------------------------------

At each time frame, an object can be described by

- `m` numbers (maximum 6), maximum three for the *location*, and
   maximum three for the *rotation*, describing the **global
   position** (location and rotation) of the object
   
- a collection of `n` numbers, describing the **local shape** of the
   object. Examples are:

	- a collection of points, describing a curve
	- a collection of points, describing a mesh
	- a collection *local locations and rotations*, describing
       assemblies of rigid bodies
	- *use your imagination!*

The full time dependent evolution of an object (a *stroke*) is
described by a curve in time **for each of the numbers above** (a
**path**). The time curves are **polynomial curves**, or **NURBS
curves**. In any case, they can be described by *coefficients* times
*basis functions*.

Short Instructions for contributions
------------------------------------------------------------

Each contribution should be done on a **feature branch** of a **forked
repository**, then a **pull request** should be issued, which will be
merged by me on mainline. Take a look at this page for a very good
explanation on possible git workflows:
https://www.atlassian.com/git/tutorials/comparing-workflows/. We will
be adopting the *forking workflow*:

- Fork the official repository (server side cloning) using the fork
  button on github
- Clone *your fork* onto your PC (it remains yours)

		git clone  https://github.com/your-user-name/ePICURE.git
		
- Add the following master remote to your clone::

		cd ePICURE
		git remote add heltai https://github.com/luca-heltai/ePICURE.git

- Make a feature branch on your PC

		git checkout heltai/master -b a_significant_name_for_your_feature
		
- commit, add tests, fix bugs, push your changes to the server,
  etc. When you push, make sure you push to *your* remote (`origin`)

		git push origin a_significant_name_for_your_feature
		
- when you think your contribution is ready (a contribution is usually
  **one single feature** with **one test**) do the following:

	- rebase interactively (`-i` switch) your changes onto the
      `master` branch of the remote `heltai`, *after* updating your
      local branches:

			# Get all changes from all remote repositories
			git fetch --all
			git rebase -i heltai/master

	- the above command will pop up an editor in which you will be
      allowed to make *squash* all your commits into a couple, more
      significant, commits. Push your changes by *forcing* the push
      (you have just changed history, so this will be allowed only if
      you force your changes)

			git push -f origin a_significant_name_for_your_feature
	
- go to github web site, select the branch you have been working on,
  and hit the button *compare and pull request*
- once I have reviewed your code, I will accept it an put it in the
  master branch


Components that should be implemented
------------------------------------------------------------

Each component should have a *well documented interface*. An example
is given for the `VectorSpace` interface, but all other interfaces
should have an equally detailed interface, which should be discussed
in the `issues` section of the git-hub repository.


**No contribution will be accepted without a small unit test!**. I'm
  using `nose` to do unit testing:
  https://nose.readthedocs.org/en/latest/. In particular it is enough
  to write a test file in the directory ./tests. If the filename
  contains the word `Test` or `test` surrounded by `_` it will be
  executed as a unit. If it contains functions whose name contains
  `test` or `Test` (surrounded by `_` or at end/begin of the name)

- An abstract python interface of type `Shape`, used to describe the
  shape of a **computable object**. [TBD: interface]
  
- An abstract python interface of type `Position`, used to describe the
  position of a **computable object**. [TBD: interface]

- An abstract python interface of type `VectorSpace`, used to describe
   *one dimensional functions* on `[0,1]`, as *coefficients* times *basis
   functions*, with access to single basis functions, their derivatives,
   their support and  the  splitting of `[0,1]` into
   sub-intervals where *each basis* is assumed to be smooth. In
   particular, the following structure should be implemented for each
   type of basis:

	- `VectorSpace.n_dofs`:  the number of basis functions that span
	  the space (*degrees of freedom*)
	- `VectorSpace.n_cells`: the number of sub-intervals of `[0,1]`
      where each basis is smooth
	- `VectorSpace.cells`: the `n_cells+1` splitting points that make
      the `cells` of the `VectorSpace`
	- `VectorSpace.basis(i)`: the ith basis function (a callable
      function)
	- `VectorSpace.basis_der(i, d)`: the d-th derivative of the i-th
      basis function (a callable function)
	- `VectorSpace.basis_span(i)`: a tuple indicating the start and
      end indices into the `cells` object where the i-th basis
      function is different from zero
	- `VectorSpace.cell_span(i)`: an array of indices containing the
      basis functions which are non zero on the i-th `cell`
	- `VectorSpace.element(c)`: a callable function,
      representing `sum(c[i] * basis[i])`, which exploits the
      *locality* of the basis functions

- The implementation of the above for

	- Power basis functions
	- Continuous Lagrange (Newton-Cotes with end points) basis
      functions, of degree `n` with `N` repetitions
	- Bernstein basis functions of degree `n` with `N` repetitions
	- B-spline basis functions with given degree and knot_vector
	- NURBS basis functions with given degree and knot_vector

- An abstract python interface of type `TimeAnalysis`, capable of
   taking derivatives and integrals of **collections of paths**, given
   an actual `VectorSpace` object, and a matrix with
   `VectorSpace.n_dofs`  rows and an arbitrary number of columns,
   representing functions  of the given `VectorSpace`

- An `Interpolation` utility, taking a `VectorSpace`, `n_dofs` points,
  and constructing a solution `c` to the problem
  `sum_over_j(c[j]*b[j](x[j])) = f(x[i])`

- A `LeastSquare` utility, doing the same but with `n>n_dofs` points
  and using least squares

- A `MassMatrix` assembler, computing the sparse matrix $M_{ij} =
  \int_0^1 v_i(x) v_j(x) dx$ exploiting locality of the `VectorSpace`
  and a `Quadrature` formula given at construction time
  
- A `StiffnessMatrix` assembler, computing the sparse stiffness matrix
  $K_{ij} = \int_0^1 v'_i(x) v'_j(x) dx$ exploiting locality of the
  `VectorSpace` and a `Quadrature` formula given at construction time

- A `ArcLengthReparametrization`, taking a matrix of coefficients
  representing a curve in a `VectorSpace`, and creating a new one via
  a `LeastSquare` fit which satisfies `c'(x) = L` for each `x`

- A `CurveFromCurvature` function, which takes a vector function
  `g(x)` representing the curvature of a curve parametrized with arc
  length, and returns the curve satisfying `c(0) = 0`, `c'(x) = L` for
  `x in [0,1]` and `c''(x) = g(x)` 

- An abstract interface `EquationOfMotion`, taking a `Shape` and
  returning a `Position`, such that some underlying equations are
  satisfied [TBD: Interface]. For this interface, we will construct

	- collection of sphere swimmers 
	- filament swimmers using resistive force theory
	- crawlers

- An abstract optimization tool, taking an `EquationOfMotion`,
  computing a `CostFunction` associated to it, imposing a
  `ConstraintEquation` on it, and optimizing the `Shape` evolution, in
  order to obtain a `Shape` that satisfies the `ConstraintEquation`
  and minimizes the `CostFunction`


- A [Blender](www.blender.org) interface

	- given two files, containing numpy nd-arrays (as created by
	   `numpy.save`), representing the time evolution of the position
	   and shape, this interface should create an animation capable of
	   visualising the time evolution of the object. We will need
	   interfaces for the following objects
	   
	    - Golestanian Axis Swimmer (`m=1, n=2`)
		- Golestanian Plane Swimmer (`m=2, n=3`)
		- Golestanian Space Swimmer (`m=6, n=4`)
		- Axis crawler (`m=1, n=as many as you wish`)
		- Plane worm (`m=2, n=as many as you wish`)
		- Space worm (`m=6, n=as many as you wish`)
		- etc.

Naming Conventions
------------------------------------------------------------

- Classes and types generally are named using uppercase letters to
  denote word beginnings (e.g. `VectorSpace`) - sometimes called
  *camel case*  -  while functions and variables use lowercase letters
  and underscores to separate words

- All file names should use lowercase letters and underscores to
  separate words

- Functions or members which return or contain the number of something
  (number of cells, degrees of freedom, etc) should start with
  `n_*`. Example: `VectorSpace.n_dofs`

- Function which set a bit or flag should start with `set_*`;
  functions which clear bits of flags should be named
  `clear_*`

- Function and variable names may not consist of only one or two
  letters, *unless* the variable is a pure counting index.

- Exceptions are used for internal parameter checking and for
  consistency checks through the `assert` mechanism. 

- Each class has to have at least 200 pages of documentation ;-)

Defensive Programming
------------------------------------------------------------

Defensive programming is a term that we use frequently when we talk
about writing code while in the mindset that errors will happen. Here,
errors can come in two ways: first, I can make a mistake myself while
writing a functions; and secondly, someone else can make a mistake
while calling my function. In either case, I would like to write my
code in such a way that errors are (i) as unlikely as possible, (ii)
that the python interpreter can already find some of the mistakes, and
(iii) that the remaining mistakes are relatively easy to find, for
example because the program aborts. Defensive programming is then a
set of strategies that make these goals more likely.

Over time, we have learned a number of techniques to this end, some of
which we list here:

- Assert preconditions on parameters: People call functions with wrong
  or nonsensical parameters, all the time. Say your function expects
  to work with two lists of the same length:

~~~ python
  
  def compare_lists(left, right):
	  # do something with left and right for each element
~~~

then make sure that you abort your program if the lists are not
compatible:

~~~ python
  
  def compare_lists(left, right):
	  assert len(left) == len(right)
	  # do something with left and right for each element
~~~

Take a look at https://wiki.python.org/moin/UsingAssertionsEffectively
for a more detailed discussion.

Notice that you can turn off assertion *automatically*: if Python is
started with the -O option, then assertions will be stripped out and
not evaluated, so you are not adding any overhead to your code. You
are only making your life easier...

License:
------------------------------------------------------------

Please see the file ./LICENSE for details


Continuous Integration Status:
------------------------

[![Build Status](https://travis-ci.org/luca-heltai/ePICURE.png)](https://travis-ci.org/luca-heltai/ePICURE)
