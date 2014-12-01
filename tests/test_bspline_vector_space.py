from interfaces.bspline_vector_space import *
from nose.tools import *

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

    # try:
    #     a.check_index(2)
    #     assert False, 'Expecting Failure!'
    # except:
    #     pass
    # assert a.basis(0)(.5) == 1
    
    # try:
    #     a.basis(0)(1.2)
    #     assert False, 'Expecting Failure!'
    # except:
    #     pass
    
