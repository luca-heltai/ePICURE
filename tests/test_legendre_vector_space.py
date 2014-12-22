import numpy as np
import matplotlib.pyplot as plt
from interfaces.legendre_vector_space import *

def test_continuous_legendre_interface():
        n = 8
        vector_space = LegendreVectorSpace(n)
        assert vector_space.n_dofs == n
        #assert vector_space.n_cells == n
        try:
            vector_space.check_index(n+2)
            assert False, 'Expecting Failure!'
        except:
            pass
        xi = vector_space.collocationJacobi(0,0,n)
        eps = 1.e-11
        for i in xrange(n):
            assert vector_space.basis(i)(xi[i]) - 1. < eps

        try:
            vector_space.basis(0)(1.2)
            assert False, 'Expecting Failure!'
        except:
            pass

        def plot_basis(i,n):
                # useful for a visual check of basis functions
                xi = np.linspace(0,1,2**8)
                f1 = plt.figure()
                for i in xrange(n):
                    plt.plot(xi,[vector_space.basis(i)(j) for j in xi],'-')
                plt.plot(vector_space.collocationJacobi(0,0,n),np.zeros(n),'ro')
                plt.show()

def test_discontinuous_legendre_interface():
        n = 8
        vector_space = DiscontinuousLegendreVectorSpace(n)
        assert vector_space.n_dofs == n
        #assert vector_space.n_cells == n
        try:
            vector_space.check_index(n+2)
            assert False, 'Expecting Failure!'
        except:
            pass
        xi = vector_space.collocationJacobi(0,0,n)
        eps = 1.e-11
        for i in xrange(n):
            assert vector_space.basis(i)(xi[i] - 1.) < eps

        try:
            vector_space.basis(0)(1.2)
            assert False, 'Expecting Failure!'
        except:
            pass
