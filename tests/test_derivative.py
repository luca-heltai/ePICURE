from interfaces import *
from utilities import *
import numpy as np
import matplotlib.pylab as plt

def check_interpolation_derivative(vs, f, df, N=20):
    print 'Checking ', type(vs).__name__

    s = np.linspace(0,1,vs.n_dofs)
    curve = f(s)

    x    = vs.element(curve)
    dx   = vs.element_der(curve,1)

    times=np.linspace(0,1,N)

    eps=1e-10
    assert (abs(x(times) - f(times))<eps).all()
    assert (abs(dx(times) - df(times))<eps).all()


def test_derivative():
    vs = UniformLagrangeVectorSpace(2)
    check_interpolation_derivative(vs, 
                                   lambda x: x, 
                                   lambda x: np.ones(np.shape(x)))
    
    vs *=3
    check_interpolation_derivative(vs, 
                                   lambda x: x,
                                   lambda x: np.ones(np.shape(x)))

    vs = UniformLagrangeVectorSpace(3)
    check_interpolation_derivative(vs, 
                                   lambda x: x**2,
                                   lambda x: 2*x)

    vs *=3
    check_interpolation_derivative(vs,
                                   lambda x: x**2,
                                   lambda x: 2*x)

    vs = UniformLagrangeVectorSpace(4)
    check_interpolation_derivative(vs,
                                   lambda x: x**3,
                                   lambda x: 3*x**2)

    vs *=3
    check_interpolation_derivative(vs,
                                   lambda x: x**3,
                                   lambda x: 3*x**2)

