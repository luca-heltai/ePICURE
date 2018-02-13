import numpy as np
from interfaces.vector_space import *
from interfaces.lagrange_vector_space import *
from utilities.matrices import *
from numpy.linalg import lstsq

def least_square_by_points(vs, points, q):
    """Returns the array c of the coefficients of the least square
    approximation of the points provided.
    """
    n = vs.n_dofs
    m = len(points)
    print (len(q), len(points), points.shape)
    assert len(q) == m, \
           'The number of points must be the same as the provided parameter subdivision'

    assert len(q) >= n, \
           'The number of points must be greater or equal to the number of dofs'
    
    M = InterpolationMatrix(vs, q)
    b = points
    c = lstsq(M, b, rcond=-1)[0]

    return c
