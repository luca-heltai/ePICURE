from interfaces.vector_space import *
import numpy as np
from numpy.polynomial.legendre import leggauss

def interpolation_matrix(vs, points):
    """Compute the Interpolation Matrix associated with the given vector space
    and the given points. This matrix can be used to solve Interpolation and

    Least Square approximations. if

    f(x) = sum(coeffs[i]*vs.basis(i)(x)), then M*coeffs = f(points)
    
    """
    # The interpolation matrix is a sparse matrix. We start by
    # constructing a full matrix, and then we transform it to sparse
    # at the end
    #
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
        M[ids,i] = vs.basis(i)(points[ids])
    return np.asmatrix(M)


def mass_matrix(vs, degree=2):
    """Create a mass matrix for the given vector space."""
    (qp, w) = leggauss(degree)
    # The default are on [-1,1]: rescale them
    qp = (qp+1)*.5
    w *= .5

    M = zeros(vs.n_dofs)
    for c in xrange(vs.n_cells):
        dofs = vs.cell_span(c)
        mloc = zeros((dofs,dofs))
        q = qp*(vs.cells[c+1]-vs.cells[c])+vs.cells[c]
        
