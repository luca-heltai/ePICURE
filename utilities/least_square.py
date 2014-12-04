import numpy as np
from interfaces.vector_space import *
from interfaces.lagrange_vector_space import *
from utilities.matrices import *
from numpy.linalg import lstsq

def least_square(f, vs, q):
    """Returns the array c of the coefficients of the least square
    approximation of the function f in the vector space vs with respect
    to the basis vs.basis.
    """
    n = vs.n_dofs
    m = len(q)
    assert m >= n, \
           'The number of points must be greater or equal to the number of dofs'
    
    M = interpolation_matrix(vs, q)
    b = f(np.array(q))
    c = lstsq(M, b)[0]

    return c
