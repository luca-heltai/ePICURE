import numpy as np

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
        q,w=((a+b)+(b-a)*qst)/2,(a+b)*wst/2
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
    

