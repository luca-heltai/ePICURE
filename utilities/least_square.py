import numpy as np
from interfaces.vector_space import *
from interfaces.lagrange_vector_space import *
from utilities.matrices import *
from numpy.linalg import solve

def least_square(f, vs, q):
    """Returns the array c of the coefficients of the least square
    approximation of the function f in the vector space vs with respect
    to the basis vs.basis.
    """
    n = vs.n_dofs
    m = len(q)
    assert m >= n, \
           'The number of points must be greater or equal to the number of dofs'
    
    B = interpolation_matrix(vs, q)
    M = B.T*B
    b = B.T*f(np.array(q).reshape((m,1)))
    c_temp = solve(M, b) #The numpy matrix type does not get along well with
    # the VectorSpace's methods!
    c = []
    for x in c_temp:
        c.append(np.asscalar(x))

    return c
