from interfaces.bspline_vector_space import *
from nose.tools import *
import numpy as np


def test_bspline_vector_space_default_constructor():
    a = BsplineVectorSpace()
    
    assert a.degree == 0
    
    assert a.knots_with_rep[0] == 0.
    assert a.knots_with_rep[1] == 1.
    assert a.knots_unique[0] == 0.
    assert a.knots_unique[1] == 1.
    
    assert a.n_knots == 2
    
    assert a.mults[0] == 1
    assert a.mults[1] == 1

    assert a.n_dofs == 1
    
    assert a.n_cells == 1

    assert a.n_dofs_per_end_point == 1

    assert a.cells[0] == 0.
    assert a.cells[1] == 1.

def test_bspline_vector_space_general_constructor():
    a = BsplineVectorSpace(2, [0.,0.,0.,.5,1.,1.,1.])
   
    assert a.degree == 2
    
    assert a.knots_with_rep[0] == 0.
    assert a.knots_with_rep[1] == 0.
    assert a.knots_with_rep[2] == 0.
    assert a.knots_with_rep[3] == 0.5
    assert a.knots_with_rep[4] == 1.
    assert a.knots_with_rep[5] == 1.
    assert a.knots_with_rep[6] == 1.
    
    assert a.knots_unique[0] == 0.
    assert a.knots_unique[1] == 0.5
    assert a.knots_unique[2] == 1.
    
    assert a.n_knots == 7

    assert a.mults[0] == 3
    assert a.mults[1] == 1
    assert a.mults[2] == 3

    assert a.n_dofs == 4

    assert a.n_cells == 2

    assert a.n_dofs_per_end_point == 1

    assert a.cells[0] == 0.
    assert a.cells[1] == 0.5
    assert a.cells[2] == 1.


def test_compute_mults():
    a = BsplineVectorSpace(3, [0.,0.,0.,0.,2.,2.,3.,5.,5.,5.,5.])
    mults = a.compute_mults(np.asarray([0.,0.,0.,0.,2.,2.,3.,5.,5.,5.,5.], np.float))

    assert len(mults) == 4
    assert mults[0] == 4
    assert mults[1] == 2
    assert mults[2] == 1
    assert mults[3] == 4


def test_compute_n_dofs():
    a = BsplineVectorSpace(3, [0.,0.,0.,0.,2.,2.,3.,5.,5.,5.,5.])
    assert a.compute_n_dofs(np.asarray([0.,2.,3.,5.], np.float),
                                       np.asarray([4,2,1,4], np.int_), 3) == 7


def test_cell_span():
    a = BsplineVectorSpace(3, [0.,0.,0.,0.,1.,2.5,5.,6.,6.,6.,6.])
    assert a.cell_span(0)[0] == 0
    assert a.cell_span(0)[1] == 1
    assert a.cell_span(0)[2] == 2
    assert a.cell_span(0)[3] == 3

    assert a.cell_span(1)[0] == 1
    assert a.cell_span(1)[1] == 2
    assert a.cell_span(1)[2] == 3
    assert a.cell_span(1)[3] == 4

    assert a.cell_span(2)[0] == 2
    assert a.cell_span(2)[1] == 3
    assert a.cell_span(2)[2] == 4
    assert a.cell_span(2)[3] == 5

    assert a.cell_span(3)[0] == 3
    assert a.cell_span(3)[1] == 4
    assert a.cell_span(3)[2] == 5
    assert a.cell_span(3)[3] == 6


def test_basis_span():
    a = BsplineVectorSpace(3, [0.,0.,0.,0.,1.,2.5,5.,6.,6.,6.,6.])
    assert a.basis_span(0) == (0.0, 1.0)
    assert a.basis_span(1) == (0.0, 2.5)
    assert a.basis_span(2) == (0.0, 5.0)
    assert a.basis_span(3) == (0.0, 6.0)
    assert a.basis_span(4) == (1.0, 6.0)
    assert a.basis_span(5) == (2.5, 6.0)
    assert a.basis_span(6) == (5.0, 6.0)


def test_find_span():
    a = BsplineVectorSpace(3, [0.,0.,0.,0.,1.,2.5,5.,6.,6.,6.,6.])
    assert a.find_span(0.) == 3
    assert a.find_span(0.1) == 3
    assert a.find_span(1.) == 3
    assert a.find_span(2.4) == 4
    assert a.find_span(2.5) == 4
    assert a.find_span(2.6) == 5
    assert a.find_span(6.) == 6

    a = BsplineVectorSpace(3, [0.,0.,0.,0.,1.,2.,2.,2.,4.,5.,6.,6.,6.,6.])
    assert a.find_span(0.) == 3
    assert a.find_span(0.1) == 3
    assert a.find_span(1.) == 3
    assert a.find_span(1.4) == 4
    assert a.find_span(2.) == 4
    assert a.find_span(4.) == 7
    assert a.find_span(5.9) == 9
    assert a.find_span(6.) == 9

def test_map_basis_cell():
    a = BsplineVectorSpace(3, [0.,0.,0.,0.,1.,2.5,5.,6.,6.,6.,6.])
    assert a.map_basis_cell(0, 3) == 0
    assert a.map_basis_cell(2, 3) == 2
    assert a.map_basis_cell(4, 5) == 2
    assert a.map_basis_cell(4, 7) == 3

    a = BsplineVectorSpace(3, [0.,0.,0.,0.,1.,2.,2.,2.,4.,5.,6.,6.,6.,6.])
    assert a.map_basis_cell(0, 3) == 0
    assert a.map_basis_cell(2, 3) == 2
    assert a.map_basis_cell(4, 7) == 0
    assert a.map_basis_cell(7, 8) == 2