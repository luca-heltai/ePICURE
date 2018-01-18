from utilities import *
from interfaces import *
import numpy as np


def check(a, b):
    print('Got shape ', b, ', expecting shape ', a)
    assert a == b, \
        'Expecting {}, got {}'.format(b, a)


def evaluate_vector_space(vs):
    """ Test the given vector space with linear coefficients. """
    print('Testing', type(vs).__name__)

    scalar_c = np.linspace(0, 1, vs.n_dofs);
    vector_c = np.arange(vs.n_dofs * 3).reshape((vs.n_dofs, 3))
    matrix_c = np.arange(vs.n_dofs * 6).reshape((vs.n_dofs, 3, 2))

    scalar_element = vs.element(scalar_c)
    vector_element = vs.element(vector_c)
    matrix_element = vs.element(matrix_c)

    # Single point
    p = .1

    # Array
    x = np.linspace(0, 1, 20)

    # Matrix
    y = np.reshape(x, (10, -1))

    # Check that shape of scalar element is the same as the input, for
    # any type.  On the other hand, the shape of a vector element,
    # should have the sum of the dimensions of the remaining
    # coefficients, plus those of the input shape.
    print('Check scalar element.')
    check(np.shape(scalar_element(p)), ())
    check(np.shape(scalar_element(x)), np.shape(x))
    check(np.shape(scalar_element(y)), np.shape(y))
    print('Check vector element.')
    check(np.shape(vector_element(p)), np.shape(vector_c[0]))
    check(np.shape(vector_element(x)), np.shape(vector_c[0]) + np.shape(x))
    check(np.shape(vector_element(y)), np.shape(vector_c[0]) + np.shape(y))
    print('Check matrix element.')
    check(np.shape(matrix_element(p)), np.shape(matrix_c[0]))
    check(np.shape(matrix_element(x)), np.shape(matrix_c[0]) + np.shape(x))
    check(np.shape(matrix_element(y)), np.shape(matrix_c[0]) + np.shape(y))

    # Same with the derivative:
    scalar_element = vs.element_der(scalar_c, 1)
    vector_element = vs.element_der(vector_c, 1)

    print('Check scalar element derivative.')
    check(np.shape(scalar_element(p)), ())
    check(np.shape(scalar_element(x)), np.shape(x))
    check(np.shape(scalar_element(y)), np.shape(y))
    print('Check vector element derivative.')
    check(np.shape(vector_element(p)), np.shape(vector_c[0]))
    check(np.shape(vector_element(x)), np.shape(vector_c[0]) + np.shape(x))
    check(np.shape(vector_element(y)), np.shape(vector_c[0]) + np.shape(y))
    print('Check matrix element.')
    check(np.shape(matrix_element(p)), np.shape(matrix_c[0]))
    check(np.shape(matrix_element(x)), np.shape(matrix_c[0]) + np.shape(x))
    check(np.shape(matrix_element(y)), np.shape(matrix_c[0]) + np.shape(y))


def test_vector_spaces():
    vs = VectorSpace()
    evaluate_vector_space(vs)
    vs = UniformLagrangeVectorSpace(2)
    evaluate_vector_space(vs)
    vs = UniformLagrangeVectorSpace(3)
    evaluate_vector_space(vs)
    vs = IteratedVectorSpace(UniformLagrangeVectorSpace(2), np.linspace(0, 1, 3))
    evaluate_vector_space(vs)
    vs = IteratedVectorSpace(UniformLagrangeVectorSpace(3), np.linspace(0, 1, 3))
    evaluate_vector_space(vs)
    vs = BsplineVectorSpace(2, [0, 0, 0, .3, .6, 1, 1, 1])
    evaluate_vector_space(vs)
