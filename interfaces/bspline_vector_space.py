import numpy as np
from interfaces.vector_space import VectorSpace


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
        
        self.degree = degree
        self.knots_with_rep = np.asarray(knots, np.float)
        self.knots_unique = np.unique(self.knots_with_rep)
        self.n_knots = self.knots_with_rep.shape[0]
        self.mults = self.compute_mults(self.knots_with_rep)

        assert ( self.n_knots - (self.degree + 1) ) > 0
        self.n_dofs = self.n_knots - (self.degree + 1)
        # self.n_dofs = self.compute_n_dofs(self.knots_unique, self.mults, self.degree)

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
        
        for i in range(1, knots_with_rep.shape[0]):
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
        for i in range(len(knots_unique)):
            n_dofs += mults[i]
        n_dofs -= (degree+1)
        return n_dofs


    def cell_span(self, i):
        """ An array of indices containing the basis functions which are non zero on the i-th cell """
        assert i >= 0
        assert i < self.n_cells

        n = 0
        for j in range(i+1):
            n += self.mults[j]
        non_zero_bases = list()
        for j in range(self.degree + 1):
            non_zero_bases.append(n - self.degree - 1 + j)

        return np.asarray(non_zero_bases, np.int_)
