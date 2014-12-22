import numpy as np
import matplotlib.pyplot as plt
from interfaces.legendre_vector_space import *

def test_continuous_legendre_interface():
        vector_space = LegendreVectorSpace(8)
        assert vector_space.n_dofs == n
        assert vector_space.n_cells == n
        try:
            vector_space.check_index(n+2)
            assert False, 'Expecting Failure!'
        except:
            pass
        xi = vector_space.collocationJacobi(0,0,n)
        for i in xrange(n):
            assert vector_space.basis(i,xi[i]) == 1

        try:
            vector_space.basis(0,1.2)
            assert False, 'Expecting Failure!'
        except:
            pass

        def plot_basis(i,n):
                # useful for a visual check of basis functions
                xi = np.linspace(0,1,2**8)
                f1 = plt.figure()
                for i in xrange(n):
                    plt.plot(xi,[vector_space.basis(i,j) for j in xi],'-')
                plt.plot(vector_space.collocationJacobi(0,0,n),np.zeros(n),'ro')
                plt.show()

test_continuous_legendre_interface()
print 'all ok'
