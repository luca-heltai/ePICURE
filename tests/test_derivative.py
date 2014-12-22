from interfaces import *
from utilities import *
import numpy as np
import matplotlib.pylab as plt

def check_linear_derivative(vs):
    print 'Checking ', type(vs).__name__

    fun=lambda t: t**2
    s = np.linspace(0,1,vs.n_dofs)
    curve = fun(s)

    x    = vs.element(curve)
    dx   = vs.element_der(curve,1)

    times=np.linspace(0,1,20)
    eps=1e-10
    plt.plot(times, dx(times))
    plt.plot(times, 2*times)
    plt.show()
             
    assert (abs(x(times) - times**2)<eps).all()
    assert (abs(dx(times) - 2*times)<eps).all()


def test_derivative():
    vs = UniformLagrangeVectorSpace(3)
    check_linear_derivative(vs)
    
    vs = IteratedVectorSpace(UniformLagrangeVectorSpace(3), np.linspace(0,1,3))
    check_linear_derivative(vs)
    
    vs = BsplineVectorSpace(2, [0,0,0,.3,.5,.8,1,1,1])
    check_linear_derivative(vs)
