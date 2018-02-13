import numpy as np
from interfaces.vector_space import *
from interfaces.lagrange_vector_space import *
from utilities.matrices import *
from numpy.linalg import solve

def interpolation(f, vs, q):
    """Returns the array c of the coefficients of the interpolation
    of the function f in the vector space vs with respect to the basis
    vs.basis, i.e. the array c such that:

    sum_{j=0}^{n}(c[j]*b[j](q[i])) = f(q[i])   for all 0<=i<=m

    where n = vs.n_dofs, m = len(q), and b[i] = vs.basis(i)
    """
    n = vs.n_dofs
    assert len(q) == n, \
           'The number of points must be equal to the number of dofs'
    
    M = interpolation_matrix(vs, q)
    b = f(np.array(q))
    c = solve(M, b)

    return c
