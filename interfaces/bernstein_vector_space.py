import numpy as np
from interfaces.vector_space import VectorSpace
from interfaces.bspline_vector_space import BsplineVectorSpace
from _Bas import basisfuns, dersbasisfuns, findspan

class BernsteinVectorSpace(BsplineVectorSpace):
    """A python interface used to describe *one dimensional Bernstein basis
    functions* on the interval given by the first and last knot.

    The base class constructs the constant vector space.
    """
    def __init__(self, degree=0):
        """ Pure interface class. It generates the constant on [0,1]. """
        assert degree >= 0
        self.degree = degree
        self.knots = np.append(np.zeros(degree+1) , np.ones(degree+1))
        self.cells = np.asarray([0,1],np.float)
        self.n_knots = len(self.knots)
        self.mults = self.compute_mults(self.knots)
        self.n_dofs = self.n_knots - (self.degree + 1)
        self.n_cells = 1
        self.n_dofs_per_end_point = 1
