import numpy as np
from interfaces.vector_space import VectorSpace
from utilities._Bas import basisfuns, dersbasisfuns

class BsplineVectorSpace(VectorSpace):
    """A python interface used to describe *one dimensional Bspline basis
    functions* onthe interval given by the first and last knot.

    The base class constructs the constant vector space.
    """
    def __init__(self, degree=0, knots=[0., 1.]):
        """ Pure interface class. It generates the constant on [0,1] if no arguments are provided. """
        assert degree >= 0
        assert len(knots) > 1
        assert knots[0] != knots[-1]

        # # Check if the knot vector is actually open.
        # if degree > 0:
        #     for i in xrange(1, degree+1):
        #         assert knots[0] == knots[i]
        #         assert knots[-1] == knots[-i]
        
        self.degree = degree
        self.knots_with_rep = np.asarray(knots, np.float)
        self.knots_unique = np.unique(self.knots_with_rep)
        self.n_knots = self.knots_with_rep.shape[0]
        self.mults = self.compute_mults(self.knots_with_rep)

        assert ( self.n_knots - (self.degree + 1) ) > 0
        #self.n_dofs = self.n_knots - (self.degree + 1)
        self.n_dofs = self.compute_n_dofs(self.knots_unique, self.mults, self.degree)

        self.n_cells = len(self.knots_unique) - 1
        self.n_dofs_per_end_point = 1
        self.cells = self.knots_unique

    def compute_mults(self, knots_with_rep):
        """Compute the multiplicity of each knot given the original knot vector.
        It returns a numpy array. It has len(knots_unique) elements.
        """
        assert len(knots_with_rep) > 1
        
        j = 1
        mults = list()
        
        for i in xrange(1, knots_with_rep.shape[0]):
            if knots_with_rep[i] == knots_with_rep[i-1]:
                j += 1
            else:
                mults.append(j)
                j = 1
        mults.append(j)

        return np.asarray(mults, np.int_)


    def compute_n_dofs(self, knots_unique, mults, degree):
        assert len(knots_unique) == len(mults)

        n_dofs = 0
        for i in xrange(len(knots_unique)):
            n_dofs += mults[i]
        n_dofs -= (degree+1)
        return n_dofs


    def cell_span(self, i):
        """ An array of indices containing the basis functions which are non zero on the i-th cell.
        They always are degree + 1."""
        assert i >= 0
        assert i < self.n_cells

        n = 0
        for j in xrange(i+1):
            n += self.mults[j]
        
        non_zero_bases = [n - self.degree - 1 + j for j in xrange(self.degree+1)]
        
        return np.asarray(non_zero_bases, np.int_)


    def basis_span(self, i):
        """Return a tuple indicating the start and end indices into the cells object where
        the i-th basis function is different from zero. Remember that a basis always spans
        degree+2 knots."""
        self.check_index(i)

        for j in xrange(len(self.cells)):
            if self.knots_with_rep[i] == self.cells[j]:
                break
        for k in xrange(len(self.cells)):
            if self.knots_with_rep[i + self.degree + 1] == self.cells[k]:
                break

        return (j, k)
        #return (self.knots_with_rep[i], self.knots_with_rep[i + self.degree + 1])

    def find_span(self, parametric_point):
        """Return the index of the knot span in which the parametric point is contained.
        The knot spans are to be intended as semi-opened (,], exept ofr the first one that 
        is closed [,]. Here knot span has to be intended with respect the complete knot 
        vector with repeted knots. The first knot has span equal to degree."""
        assert parametric_point >= self.knots_unique[0]
        assert parametric_point <= self.knots_unique[-1]

        i = self.mults[0]-1
        while parametric_point > self.knots_with_rep[i+1]:
            i += 1
        return i


    def map_basis_cell(self, i, knot_interval):
        """This method returns the index of the cell_span vector that corresponds to the 
        i-th basis function. The cell_span takes the cell with respect the knots_unique, 
        while the knot_interval is with respect the knots_with_rep. So we have to sum 
        the multiplicity of the knots -1 untill the first knot of the current interval."""
        knot = self.knots_with_rep[knot_interval]
        n = 0
        for j in xrange(len(self.knots_unique)):
            if self.knots_unique[j] == knot:
                n = j
        
        summation = 0
        for j in xrange(n+1):
            summation += self.mults[j]-1

        non_zero_bases = self.cell_span(knot_interval - summation)
        for j in xrange(self.degree+1):
            if non_zero_bases[j] == i:
                return j


    def basis(self, i):
        """The ith basis function (a callable function)"""
        self.check_index(i)
        t = self.basis_span(i)
        # If the point is outside the support of the i-th basis function it returns 0.
        g = lambda x: basisfuns(self.find_span(x), self.degree, x, 
            self.knots_with_rep)[self.map_basis_cell(i, self.find_span(x))] if x >= self.cells[t[0]] and x <= self.cells[t[1]] else 0
        return np.frompyfunc(g, 1, 1)


    def basis_der(self, i, d):
        """The d-th derivative of the i-th basis function (a callable function)"""
        self.check_index(i)
        t = self.basis_span(i)
        # If the point is outside the support of the i-th basis function it returns 0.
        # We take the d-th row of the matrix returned by dersbasisfuns
        g = lambda x: dersbasisfuns(self.degree, self.knots_with_rep, x, self.find_span(x),
            d)[d][self.map_basis_cell(i, self.find_span(x))] if x >= self.cells[t[0]] and x <= self.cells[t[1]] else 0
        return np.frompyfunc(g, 1, 1)



