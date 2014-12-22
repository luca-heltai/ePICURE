import numpy as np
import scipy.special
from interfaces.lagrange_vector_space import LagrangeVectorSpace


class LegendreVectorSpace(LagrangeVectorSpace):
    """A python interface used to describe *one dimensional lagrange basis functions
    on n legendre collocation points*.
    """
    def __init__(self, n, continuous=True):
         LagrangeVectorSpace.__init__(self, self.collocationJacobi(0,0,n), continuous)

    def collocationJacobi(self, alpha, beta, n):
         Jroots = np.zeros(n)
         # Jacobi polynomials are defined on (-1,1): rescale to (0,1)
         assert (alpha == 0 and beta == 0), 'Using non-Legendre Jacobi polynomials. \
                         Careful: weak formulation needs weights.'
         assert (n > 2), 'Polynomial degree too low for Legendre polynomials.\
                         Use linear Lagrange instead.'
         Jroots[1:-1] = 0.5*(scipy.special.j_roots(n-2, alpha, beta)[0]+1)
         Jroots[-1] = 1.
         return Jroots
