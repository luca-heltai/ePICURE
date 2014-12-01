from interfaces.bspline_vector_space import *
from nose.tools import *
import numpy as np

def test_bspline_vector_space_default_constructor():
    a = BsplineVectorSpace()
    
    assert a.degree == 0
    
    assert a.knots_with_rep[0] == 0.
    assert a.knots_with_rep[1] == 1.
    assert a.knots_unique[0] == 0.
    assert a.knots_unique[0] == 0.
    
    assert a.n_knots == 2
    
    assert a.mults[0] == 1
    assert a.mults[1] == 1

    assert a.n_dofs == 1
    

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

