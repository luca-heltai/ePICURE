from interfaces import *
from utilities import *
import numpy as np
import matplotlib.pylab as plt

def test_massmatrix():
    for i in range(10):
        for j in range(2,i):
            vs=IteratedVectorSpace(UniformLagrangeVectorSpace(j), np.linspace(0,1,i))
            v=np.ones(vs.n_dofs)
            for k in range(2,10):
                t=massmatrix(vs,k)
                print v.dot(t.dot(v))

test_massmatrix()
