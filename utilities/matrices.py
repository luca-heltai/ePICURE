from interfaces.vector_space import *
import numpy as np
from numpy.polynomial.legendre import leggauss

def interpolation_matrix(vs, points, d=0):
    """Compute the Interpolation Matrix associated with the given vector space
    and the given points. This matrix can be used to solve Interpolation and
    Least Square approximations. The matrix which is returned is defined as

    f(x) = sum(coeffs[i]*vs.basis_der(i,d)(x)), then M*coeffs = f(points)
    
    """
    # M_ij := v_j(q_i)
    col=points.reshape((-1,1))
    M = np.zeros((len(points), vs.n_dofs), order='F')
    for i in xrange(vs.n_dofs):
        # advanced slicing... only compute those points which are
        # different from zero
        ia,ib = vs.basis_span(i)
        a = vs.cells[ia]
        b = vs.cells[ib]
        ids = (a<=points) & (points<=b)
        if d == 0:
            M[ids,i] = vs.basis(i)(points[ids])
        else:
            M[ids,i] = vs.basis_der(i,d)(points[ids])
    return np.asmatrix(M)


        
