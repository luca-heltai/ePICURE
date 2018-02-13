from interfaces import *
from utilities import *
import numpy as np


def test_mass_matrix():
    print ("Mass matrix test")
    print ("This should take ~30s in a core i3")
    for i in range(10):
        for j in range(2,i):
            vs = IteratedVectorSpace(UniformLagrangeVectorSpace(j), np.linspace(0,1,num=i))
            v = np.ones(vs.n_dofs)
            for k in range(2, 10):
                t = MassMatrix(vs, k)
                assert abs(v.dot(t.dot(v))-1.0) < 10**(-10)
                csc = MassMatrix(vs, k, matrix_format="CSC")
                csr = MassMatrix(vs, k, matrix_format="CSR")
                assert (t == csc).all() and (t == csr).all()


def test_stiffness_matrix():
    print ("Stiffness matrix test")
    print("This should take ~30s in a core i3")
    for i in range(10):
        for j in range(2, i):
            vs = IteratedVectorSpace(UniformLagrangeVectorSpace(j), np.linspace(0, 1, num=i))
            v = np.ones(vs.n_dofs)
            for k in range(2, 10):
                t = StiffnessMatrix(vs, k)
                csc = StiffnessMatrix(vs, k, matrix_format="CSC")
                csr = StiffnessMatrix(vs, k, matrix_format="CSR")
                assert (t == csc).all() and (t == csr).all()

#test_mass_matrix()
