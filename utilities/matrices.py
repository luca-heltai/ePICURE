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

def massmatrix(vs,nq,quadfun=np.polynomial.legendre.leggauss):
    """A MassMatrix assembler, computing the sparse matrix 
    
    $M_{ij} =\int_0^1 v_i(x) v_j(x) dx$
    vs: A vector_space type opject.
    nq: Number of quadrature points used for each individual integration.
    quadfun: a function that receives a number and returns a touple of 
                numpy arrays with the values of quadrature points and 
                the weights in the interval [-1;1]
    """
    (qst,wst)=quadfun(nq)
    cum=np.zeros((vs.n_dofs,vs.n_dofs))
    for lcell in xrange(vs.n_cells):
        a,b=vs.cells[lcell],vs.cells[lcell+1]
        q,w=((a+b)+(b-a)*qst)/2,(b-a)*wst/2
        fnts=vs.cell_span(lcell)
        nfnts=len(fnts)
        fq=np.zeros((nfnts,nq))
        for i in range(nfnts):                                          #funtion evaluation is a heavy operation do it one time
            fq[i]=vs.basis(fnts[i])(q)
        out=np.dot(fq*w,fq.transpose())                                 #all in only one step
        for rout, rcum in zip(out, fnts):
            for vout, ccum in zip(rout, fnts):
                cum[rcum,ccum]+=vout
    return cum
