import numpy as np
from interfaces.vector_space import VectorSpace
from igakit.igalib import bsp
from igakit.cad import line


class BsplineVectorSpace(VectorSpace):
    """A python interface used to describe *one dimensional Bspline basis
    functions* onthe interval given by the first and last knot.

    The base class constructs the constant vector space.
    """
    def __init__(self, degree=0, knots=[]):
        """ Pure interface class. It generates the constant on [0,1] if no arguments are provided. """
        assert degree >= 0

        self.degree = degree
        if len(knots) == 0:
            self.knots = (self.degree + 1) * [0] + (self.degree + 1) * [1]
        else:
            assert knots[0] != knots[-1]
            self.knots = knots
        self.cells = np.unique(self.knots)
        self.n_knots = len(self.knots)
        self.mults = self.compute_mults(self.knots)

        assert ( self.n_knots - (self.degree + 1) ) > 0
        self.n_dofs = self.n_knots - (self.degree + 1)

        self.n_cells = len(self.cells) - 1
        self.n_dofs_per_end_point = 1

    def compute_mults(self, knots):
        """Compute the multiplicity of each cell boundary, given the original
        knot vector.  It returns a numpy array with the same length as cells.
        """
        assert len(knots) > 1
        
        j = 1
        mults = list()
        
        for i in range(1, len(knots)):
            if knots[i] == knots[i-1]:
                j += 1
            else:
                mults.append(j)
                j = 1
        mults.append(j)

        return np.asarray(mults, np.int_)


    def cell_span(self, i):
        """An array of indices containing the basis functions which are non zero on
        the i-th cell.  If the knot span is closed, there are always
        degree+1 non zero basis functions. In other cases, the first and last cell
        may have a different number. 
        """
        list_of_dofs = self.internal_cell_span(i)
        return  [ i for i in list_of_dofs if i in range(0, self.n_dofs) ]


    def internal_cell_span(self, i):
        """ An array of indices containing the basis functions which are non zero on the i-th cell.
        They always are degree + 1."""
        assert i >= 0
        assert i < self.n_cells

        n = 0
        for j in range(i+1):
            n += self.mults[j]
        
        non_zero_bases = [n - self.degree - 1 + j for j in range(self.degree+1)]
        
        return np.asarray(non_zero_bases, np.int_)

    def basis_span(self, i):
        """Return a tuple indicating the start and end indices into the cells object where
        the i-th basis function is different from zero. Remember that a basis always spans
        degree+2 knots."""
        self.check_index(i)
        j = np.where(self.cells == self.knots[i])[0][0]
        k = np.where(self.cells == self.knots[i+1+self.degree])[0][0]
        return (j, k)


    def find_span(self, parametric_point):
        """Return the index of the knot span in which the parametric point is contained.
        The knot spans are to be intended as semi-opened [,), except on the last one that 
        is closed [,]. Here knot span has to be intended with respect the complete knot 
        vector with repeted knots. """
        assert parametric_point >= self.cells[0]
        assert parametric_point <= self.cells[-1]
        if parametric_point == self.cells[0]:
            return (self.mults[0]-1)
        elif parametric_point == self.cells[-1]:
            return (self.n_knots-self.mults[-1]-1)
        else:
            return (np.where(np.asarray(self.knots) <= parametric_point)[0][-1])


    def map_basis_cell(self, i, knot_interval):
        """This method returns the index of the cell_span vector that corresponds to the 
        i-th basis function. The cell_span takes the cell with respect the cells, 
        while the knot_interval is with respect the knots. So we have to sum 
        the multiplicity of the knots -1 untill the first knot of the current interval."""
        knot = self.knots[knot_interval]
        n = 0
        for j in range(len(self.cells)):
            if self.cells[j] == knot:
                n = j
        
        summation = 0
        for j in range(n+1):
            summation += self.mults[j]-1

        non_zero_bases = self.internal_cell_span(knot_interval - summation)
        for j in range(self.degree+1):
            if non_zero_bases[j] == i:
                return j
        # By default, return degree
        return self.degree


    def basis(self, i):
        """The ith basis function (a callable function)"""
        self.check_index(i)
        t = self.basis_span(i)
        # If the point is outside the support of the i-th basis function it returns 0.
        #g = lambda x: basisfuns(self.find_span(x), self.degree, x,
        #    self.knots)[self.map_basis_cell(i, self.find_span(x))] if x >= self.cells[t[0]] and x <= self.cells[t[1]] else 0
        g = lambda x: bsp.EvalBasisFuns(self.degree, self.knots, x)[self.map_basis_cell(i, self.find_span(x))]\
            if x >= self.cells[t[0]] and x <= self.cells[t[1]] else 0
        return np.frompyfunc(g, 1, 1)


    def basis_der(self, i, d):
        """The d-th derivative of the i-th basis function (a callable function)"""
        self.check_index(i)
        t = self.basis_span(i)
        # If the point is outside the support of the i-th basis function it returns 0.
        # We take the d-th row of the matrix returned by dersbasisfuns
        #g = lambda x: dersbasisfuns(self.degree, self.knots, x, self.find_span(x),
        #    d)[d][self.map_basis_cell(i, self.find_span(x))] if x >= self.cells[t[0]] and x <= self.cells[t[1]] else 0
        assert d <= 3, "igakit supports up to the third derivative"
        g = lambda x: bsp.EvalBasisFunsDers(self.degree, self.knots, x, d)[d][self.map_basis_cell(i, self.find_span(x))]\
            if x >= self.cells[t[0]] and x <= self.cells[t[1]] else 0
        return np.frompyfunc(g, 1, 1)


    def greville(self):
        """Return the set of greville abscissa associated with this Bspline"""
        k = self.knots
        g = np.zeros(self.n_dofs)
        for i in range(self.n_dofs):
            g[i] = np.average(k[i+1:i+self.degree+1])
        return g
    
class UniformBsplineVectorSpace(BsplineVectorSpace):
    def __init__(self, p, nel):
        nbasis = max(2, nel+1)
        knots = np.r_[(p)*[0], np.linspace(0,1, nbasis), (p)*[1]] 
        BsplineVectorSpace.__init__(self, p, knots)

