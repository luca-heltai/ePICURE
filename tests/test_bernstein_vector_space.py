from interfaces.bernstein_vector_space import *
from nose.tools import *
import numpy as np


def test_bernstein_vector_space_default_constructor():
    b = BernsteinVectorSpace()

    assert b.degree == 0

    assert b.knots[0] == 0.
    assert b.knots[1] == 1.

    assert b.cells[0] == 0.
    assert b.cells[1] == 1.

    assert b.n_knots == 2

    assert b.mults[0] == 1
    assert b.mults[1] == 1

    assert b.n_dofs == 1

    assert b.n_cells == 1

    assert b.n_dofs_per_end_point == 1


def test_bernstein_vector_space_general_constructor():
    b = BernsteinVectorSpace(2)

    assert b.degree == 2

    assert b.knots[0] == 0.
    assert b.knots[1] == 0.
    assert b.knots[2] == 0.
    assert b.knots[3] == 1.
    assert b.knots[4] == 1.
    assert b.knots[5] == 1.

    assert b.cells[0] == 0.
    assert b.cells[1] == 1.

    assert b.n_knots == 6

    assert b.mults[0] == 3
    assert b.mults[1] == 3

    assert b.n_dofs == 3

    assert b.n_cells == 1

    assert b.n_dofs_per_end_point == 1


def test_compute_mults():
    b = BernsteinVectorSpace(3)
    mults = b.compute_mults(b.knots)
    assert len(mults) == 2
    assert mults[0] == 4
    assert mults[1] == 4


def test_cell_span():
    b = BernsteinVectorSpace(3)
    assert b.cell_span(0)[0] == 0
    assert b.cell_span(0)[1] == 1
    assert b.cell_span(0)[2] == 2
    assert b.cell_span(0)[3] == 3


def test_basis():
    b = BernsteinVectorSpace(3)
    assert b.basis(0)(0) == 1.
    assert b.basis(0)(1) == 0.


def test_basis_der():
    b = BernsteinVectorSpace(3)
    assert b.basis_der(0,0)(0) == 1.
    assert b.basis_der(0,1)(0) == -3.0
    assert b.basis_der(0,2)(0) == 6.0
    assert b.basis_der(0,3)(0) == -6.0
    #assert b.basis_der(0,4)(0) == 0.

