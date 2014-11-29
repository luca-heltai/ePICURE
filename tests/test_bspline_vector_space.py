from interfaces.bspline_vector_space import *
from nose.tools import *

def test_bspline_vector_space_default_constructor():
    a = BsplineVectorSpace()
    assert a.degree == 0
    assert a.knots[0] == 0.
    assert a.knots[1] == 1.
    assert a.knots_unique[0] == 0.
    assert a.knots_unique[0] == 0.
    assert a.n_knots == 2
    

def test_bspline_vector_space_general_constructor():
    a = BsplineVectorSpace(2, [0.,0.,0.,.5,1.,1.,1.])
    assert a.degree == 2
    assert a.knots[0] == 0.
    assert a.knots[1] == 0.
    assert a.knots[2] == 0.
    assert a.knots[3] == 0.5
    assert a.knots[4] == 1.
    assert a.knots[5] == 1.
    assert a.knots[6] == 1.
    
    assert a.knots_unique[0] == 0.
    assert a.knots_unique[1] == 0.5
    assert a.knots_unique[2] == 1.
    assert a.n_knots == 7

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
    
