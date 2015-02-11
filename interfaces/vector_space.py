import numpy as np
from itertools import product as cartesian_product

class VectorSpace(object):
    """An abstract python interface used to describe *one dimensional
    functions* on `[a,b]`, as *coefficients* times *basis functions*,
    with access to single basis functions, their derivatives, their
    support and the splitting of `[a,b]` into sub-intervals where
    *each basis* is assumed to be smooth.

    This class provides an abstract interface to describe VectorSpaces
    of functions from `[a,b]` into general `R^N = R^mxnxd...`.

    The main attribute of this class is the `element` method, which
    takes as input a numpy array of length VectorSpace.n_dofs

    len(C) == VectorSpace.n_dofs # True    
    g = VectorSpace.element(C)
    
    In general, this results in a callable function g, such that

    shape( g(x) ) = shape(C[0])+shape(x)
    
    g := C[i,j,k, ...]*VectorSpace.basis(i) # sum over i is implicit

    If you evaluate g at a nparray x, then
    
    C = g(x)
    C[j,k,...,l,m,...] := C[i,j,k, ...]*VectorSpace.basis(i)(x[l,m, ...])
    
    where again a sum is implicit in i.
    """
    def __init__(self, a=0.0, b=1.0):
        """ Pure interface class. It generates the constant on [0,1]. """
        self.n_dofs = 1
        self.n_dofs_per_end_point = 0
        self.n_cells = 1
        self.cells = np.array([a, b])
        self.degree = 0

    def check_index(self, i):
        assert i< self.n_dofs, \
            'Trying to access base {}, but we only have {}'.format(i, self.n_dofs)

    def basis(self, i):
        self.check_index(i)
        return lambda x: np.ones(np.shape(x))

    def basis_der(self, i, d):
        self.check_index(i)
        return lambda x: np.zeros(np.shape(x))

    def basis_span(self, i):
        """ VectorSpace.basis_span(i): a tuple indicating the start and end indices into the cells object where the i-th basis function is different from zero"""
        self.check_index(i)
        return (0, 1)

    def cell_span(self, i):
        """ An array of indices containing the basis functions which are non zero on the i-th cell """
        self.check_index(i)
        return [0]

    def eval(self, c, x):
        assert len(c) == self.n_dofs, \
            'Incompatible vector. It should have length {}. It has lenght {}'.format(self.n_dofs, len(c))
        # Find the cell where x lies:
        sy = np.shape(c[0])+np.shape(x)
        y = np.zeros(sy)
        # TBD: make this function aware of locality
        for i in xrange(self.n_dofs):
            y += np.outer(c[i], self.basis(i)(x)).reshape(sy)
        return y
    
    def eval_der(self, c, d, x):
        assert len(c) == self.n_dofs, \
          'Incompatible vector. It should have length {}. It has lenght {}'.format(self.n_dofs, len(c))
        # Find the cell where x lies:
        sy = np.shape(c[0])+np.shape(x)
        y = np.zeros(sy)
        # TBD: make this function aware of locality
        for i in xrange(self.n_dofs):
            y += np.outer(c[i], self.basis_der(i,d)(x)).reshape(sy)
        return y
        
    def element(self, c):
        """VectorSpace.element(c): a callable function, representing sum(c[i] *
        basis[i]), which exploits the locality of the basis functions
        """
        return lambda x: self.eval(c, x)
    
    def element_der(self, c, d):
        """VectorSpace.element(c): a callable function, representing sum(c[i] *
        basis[i]), which exploits the locality of the basis functions
        """
        return lambda x: self.eval_der(c, d, x)

    def print_info(self):
        print '============================================================'
        print 'Name: '+type(self).__name__
        print 'N dofs: {}, N cells: {}, \nCell boundaries: {}'.format(self.n_dofs, self.n_cells, self.cells)
        print 'Shared dofs on cell boundaries: {}'.format(self.n_dofs_per_end_point)
        for i in xrange(self.n_cells):
            print '------------------------------------------------------------'
            print 'Cell {}: [{},{}]'.format(i, self.cells[i], self.cells[i+1])
            print 'Nonzero basis: {}'.format(self.cell_span(i))
        print '------------------------------------------------------------'

    def __imul__(self, n):
        """This allows the construction of repeated VectorSpaces by simply
        multiplying the VectorSpace with an integer. """
        return RepeatedVectorSpace(self, n)
        

class TensorProduct(VectorSpace):
    def __init__(self, space_list):
        self.space_list = space_list
        shape = []
        self.n_dofs = 1
        self.n_cells = 1
        for i in space_list:
            self.n_dofs *= i.n_dofs
            self.n_cells *= i.n_cells
            shape.append(i.n_dofs)
        self.shape = tuple(shape)
        
    def basis(self, i):
        self.check_index(i)
        indx = np.unravel_index(i, self.shape)
        def basis_i(x):
            output = 1
            for j in range(len(self.shape)):
                output *= self.space_list[j].basis(indx[j])(x[j])
            return output
        return basis_i

    def basis_der(self, i, d):
        self.check_index(i)
        indx = np.unravel_index(i, self.shape)
        def basis_der_i(x):
            output = 1
            for j in range(len(self.shape)):
                output *= self.space_list[j].basis_der(indx[j],d)(x[j])
            return output
        return basis_der_i

    def basis_span(self, i):
        self.check_index(i)
        indx = np.unravel_index(i, self.shape)
        min_indices = []
        max_indices = []
        for j in xrange(len(self.shape)):
            min_indices.append(np.min(self.space_list[j].basis_span(indx[j])))
            max_indices.append(np.max(self.space_list[j].basis_span(indx[j])))
        i1 = np.ravel_multi_index(min_indices, self.shape)
        i2 = np.ravel_multi_index(max_indices, self.shape)
        return (i1, i2)

    def cell_span(self, i):
        cell_span_list = []
        indx = np.unravel_index(i, self.shape)
        output = []
        cell_for_axis = [ self.space_list[j].cell_span(indx[j]) for j in range(len(indx))]
        for multi_indx in cartesian_product(*cell_for_axis):
            output.append(np.ravel_multi_index(multi_indx, self.shape))
        return output




class AffineVectorSpace(VectorSpace):
    """Affine transformation of a VectorSpace. Given a vector space, returns
    its affine transformation between a and b
    """
    def __init__(self, vs, a=0.0, b=1.0):
        """ Scale the vector space from its original span to [a,b]"""
        self.vs = vs
        self.n_dofs = vs.n_dofs
        self.n_dofs_per_end_point = vs.n_dofs_per_end_point
        self.n_cells = vs.n_cells
        a0 = vs.cells[0]
        b0 = vs.cells[-1]
        self.J = (b-a)/(b0-a0)
        self.cells = (vs.cells-a0)/(b0-a0)*(b-a) + a
        self.a0 = a0
        self.b0 = b0
        self.a = a
        self.b = b
        self.degree = vs.degree

    def reset(self, a=0.0, b=1.0):
        """Make the affine transformation on the new [a,b] interval."""
        a0 = self.a0
        b0 = self.b0
        self.a = a
        self.b = b
        self.J = (b-a)/(b0-a0)
        self.cells = (self.vs.cells-a0)*self.J + a
        
    def pull_back(self, x):
        """Transform x from [a,b] to [a0, b0]"""
        return (x-self.a)/self.J+self.a0
        
    def push_forward(self, x0):
        """Transform x from [a0,b0] to [a, b]"""
        return (x0-self.a0)*self.J+self.a
    
    def basis(self, i):
        return lambda x: self.vs.basis(i)(self.pull_back(x))

    def basis_der(self, i, d):
        return lambda x: self.vs.basis_der(i,d)(self.pull_back(x))/(self.J**d)
 
    def basis_span(self, i):
        """VectorSpace.basis_span(i): a tuple indicating the start and end indices
        into the cells object where the i-th basis function is different from
        zero
        """
        return self.vs.basis_span(i)

    def cell_span(self, i):
        """ An array of indices containing the basis functions which are non zero on the i-th cell """
        return self.vs.cell_span(i)



class IteratedVectorSpace(VectorSpace):
    """Construct an iterated version of the given original VectorSpace
    """
    def __init__(self, vsext, span):
        """ Construct this vectorspace as a replica of vs in each of the span interval."""
        assert len(span) > 1, \
          'Expecting span to have at least two entries! Found {}'.format(len(span))
          
        self.vs = AffineVectorSpace(vsext, 0, 1)
        vs = self.vs
        self.span = span
        # Start by assuming discontinuous functions
        self.n_dofs = vs.n_dofs*(len(span)-1)
        self.n_dofs_per_end_point = vs.n_dofs_per_end_point
        if self.n_dofs_per_end_point > 0:
            self.n_dofs -= (len(span)-2)*vs.n_dofs_per_end_point
        self.n_cells = vs.n_cells*(len(span)-1)

        # All cells have already been set to 0,1. Remove the last
        # element, which is repeated
        loc_cells = vs.cells[0:-1]
        self.cells = np.array([])
        for i in xrange(len(span)-1):
            self.cells = np.append(self.cells, loc_cells*(span[i+1]-span[i])+span[i])
        self.cells = np.append(self.cells, [span[-1]])
        assert len(self.cells)-1 == self.n_cells, \
          "Internal error! {} != {}".format(len(self.cells)-1, self.n_cells)
        self.degree = vs.degree

    def local_to_global(self, base, i):
        """Return global index of ith local dof on cell c"""
        self.vs.check_index(i)
        # Keep track of the fact that the last n_dofs_per_end_point
        # are shared between consecutive basis
        return self.vs.n_dofs*base+i-base*self.n_dofs_per_end_point
        

    def global_to_local(self, i):
        """Given a global index, return the local base index (or indices), and the
        local dof (or dofs) to which this degree of freedom refers to as
        pairs of [ (local_vspace_index, local_index),
        (local_vspace_index, local_index)].  This has nothing to do with
        the cells indices. Those are returned by basis_span.
        """
        assert self.vs.n_dofs-self.n_dofs_per_end_point > 0, \
          'Internal error! {} ! > 0'.format(self.vs.n_dofs-self.n_dofs_per_end_point) 

        n_unique = (self.vs.n_dofs-self.n_dofs_per_end_point)
        b = int(np.floor(i/n_unique))
        j = np.mod(i, n_unique)
        ret = []
        if j < self.n_dofs_per_end_point and b>0:
            ret += [(b-1, self.vs.n_dofs-(self.n_dofs_per_end_point-j))]
        if b< len(self.span)-1:
            ret += [(b,j)]
        return ret

    def eval_basis(self, i, xx):
        self.check_index(i)
        pairs = self.global_to_local(i)
        x = np.squeeze(xx)
        y = np.squeeze(np.array(x*0))
        span = self.span
        for p in pairs:
            a = span[p[0]]
            b = span[p[0]+1]
            vs = self.vs
            vs.reset(a, b)
            if b == self.cells[-1]:
                b += 1
            ids = np.array( (a<=x) & (x<b) )
            y[ids] = vs.basis(p[1])(x[ids])
        return y
            
            
    def eval_basis_der(self, i, d, xx):
        self.check_index(i)
        pairs = self.global_to_local(i)
        x = np.squeeze(xx)
        y = np.squeeze(np.array(x*0))
        span = self.span
        for p in pairs:
            a = span[p[0]]
            b = span[p[0]+1]
            vs = self.vs
            vs.reset(a, b)
            if b == self.cells[-1]:
                b += 1
            ids = np.array( (a<=x) & (x<b) )
            y[ids] = vs.basis_der(p[1], d)(x[ids])
        return y
    
    def basis(self, i):
        return lambda x: self.eval_basis(i,x)

    def basis_der(self, i, d):
        return lambda x: self.eval_basis_der(i,d,x)

    def basis_span(self, i):
        """VectorSpace.basis_span(i): a tuple indicating the start and end indices
        into the cells object where the i-th basis function is different from
        zero
        """
        self.check_index(i)
        nbasis = len(self.span)-1
        pairs = self.global_to_local(i)
        start = self.vs.n_cells*pairs[0][0]+self.vs.basis_span(pairs[0][1])[0]
        end   = self.vs.n_cells*pairs[-1][0]+self.vs.basis_span(pairs[-1][1])[1]
        return (start, end)

    def cell_span(self, i):
        """ An array of indices containing the basis functions which are non zero on the i-th cell """
        b = i/self.vs.n_cells
        j = np.mod(i, self.vs.n_cells)
        startid = b*(self.vs.n_dofs-self.vs.n_dofs_per_end_point)
        local_ids = self.vs.cell_span(j)
        return [id+startid for id in local_ids]

class RepeatedVectorSpace(IteratedVectorSpace):
    """Construct an IteratedVectorSpace with uniform repetitions."""
    def __init__(self, vs, n):
        IteratedVectorSpace.__init__(self, vs, np.linspace(0,1,n+1))

