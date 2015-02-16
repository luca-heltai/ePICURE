from interfaces.vector_space import *
from interfaces.lagrange_vector_space import *
from numpy import array

def test_tensor_product_1Element():
    a = VectorSpace()
    ts = TensorProduct([a])
    assert ts.n_dofs == 1
    assert ts.n_cells == 1
    try:
        ts.check_index(2)
        assert False, 'Expecting Failure!'
    except:
        pass
    print ts.basis(0)(array([.5]))
    assert ts.basis(0)(array([.5])) == 1
    
    try:
        ts.basis(0)(1.2)
        assert False, 'Expecting Failure!'
    except:
        pass

def test_tensor_product_2Elements():
    a = VectorSpace()
    ts = TensorProduct([a,a])
    assert ts.n_dofs == 1
    assert ts.n_cells == 1
    assert ts.basis(0)(array([.5,.5])) == 1

def test_tensor_product_3Elements():
    a = VectorSpace()
    vs = UniformLagrangeVectorSpace(7)
    ts = TensorProduct([a,vs,a])
    assert ts.n_dofs == vs.n_dofs
    assert ts.n_cells == vs.n_cells

def test_tensor_product_evaluations():
    vs1 = UniformLagrangeVectorSpace(2)
    vs2 = UniformLagrangeVectorSpace(3)
    ts = TensorProduct([vs1,vs2])
    print ts.n_dofs
    print vs1.n_dofs
    print vs2.n_dofs
    assert ts.n_dofs == vs1.n_dofs * vs2.n_dofs
    assert ts.n_cells == vs1.n_cells * vs2.n_cells

    test_point = array([.5, .5])

    assert ts.basis(0)(test_point) == vs1.basis(0)(.5) * vs2.basis(0)(.5)
    assert ts.basis(1)(test_point) == vs1.basis(0)(.5) * vs2.basis(1)(.5)
    assert ts.basis(2)(test_point) == vs1.basis(0)(.5) * vs2.basis(2)(.5)
    assert ts.basis(3)(test_point) == vs1.basis(1)(.5) * vs2.basis(0)(.5)
    assert ts.basis(4)(test_point) == vs1.basis(1)(.5) * vs2.basis(1)(.5)
    assert ts.basis(5)(test_point) == vs1.basis(1)(.5) * vs2.basis(2)(.5)

    assert ts.basis_der(0,1)(test_point) == vs1.basis_der(0,1)(.5) * vs2.basis_der(0,1)(.5)
    assert ts.basis_der(1,1)(test_point) == vs1.basis_der(0,1)(.5) * vs2.basis_der(1,1)(.5)
    assert ts.basis_der(2,1)(test_point) == vs1.basis_der(0,1)(.5) * vs2.basis_der(2,1)(.5)
    assert ts.basis_der(3,1)(test_point) == vs1.basis_der(1,1)(.5) * vs2.basis_der(0,1)(.5)
    assert ts.basis_der(4,1)(test_point) == vs1.basis_der(1,1)(.5) * vs2.basis_der(1,1)(.5)
    assert ts.basis_der(5,1)(test_point) == vs1.basis_der(1,1)(.5) * vs2.basis_der(2,1)(.5)

def test_tensor_product_basis_span():
    vs1 = IteratedVectorSpace(UniformLagrangeVectorSpace(3), np.linspace(0,1,3))
    vs2 = IteratedVectorSpace(UniformLagrangeVectorSpace(2), np.linspace(0,1,4))
    ts = TensorProduct([vs1,vs2])
    a = vs1.basis_span(1)[0]
    b = vs2.basis_span(3)[0]
    c = a * vs2.n_dofs + b
    assert ts.basis_span(1*vs2.n_dofs + 3)[0] == c
    a = vs1.basis_span(1)[1]
    b = vs2.basis_span(3)[1]
    c = a * vs2.n_dofs + b
    assert ts.basis_span(1*vs2.n_dofs + 3)[1] == c


def test_tensor_product_cell_span():
    vs1 = IteratedVectorSpace(UniformLagrangeVectorSpace(2), np.linspace(0,1,3))
    vs2 = IteratedVectorSpace(UniformLagrangeVectorSpace(4), np.linspace(9,11,4))
    ts = TensorProduct([vs1,vs2])
    
    assert len(ts.cell_span(1*vs2.n_dofs + 2)) == len(vs1.cell_span(1)) * len(vs2.cell_span(2))
