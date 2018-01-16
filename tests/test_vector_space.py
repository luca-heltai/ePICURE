from interfaces.vector_space import *

def test_vector_space_interface():
    a = VectorSpace()
    assert a.n_dofs == 1
    assert a.n_cells == 1
    try:
        a.check_index(2)
        assert False, 'Expecting Failure!'
    except:
        pass
    assert a.basis(0)(.5) == 1
    
    try:
        a.basis(0)(1.2)
        assert False, 'Expecting Failure!'
    except:
        pass
