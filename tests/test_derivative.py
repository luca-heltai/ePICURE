from interfaces import *
from utilities import *
import numpy as np
import matplotlib.pylab as plt

def check_interpolation_derivative(vs, f, df, d=1, N=20):
    print 'Checking ', type(vs).__name__ , 'derivative', d

    s = np.linspace(vs.cells[0], vs.cells[-1], vs.n_dofs)
    curve = f(s)

    x    = vs.element(curve)
    dx   = vs.element_der(curve,d)

    times=np.linspace(vs.cells[0], vs.cells[-1], N)

    eps=1e-10
    c = vs.cells[-1]
    print 'Df^',d,'(', c, ') =', df(c), ' - approx - ', dx(c)
    assert (abs(x(times) - f(times))<eps).all()
    assert (abs(dx(times) - df(times))<eps).all()


def test_derivative():
    f =  [lambda x: x,
          lambda x: x**2,
          lambda x: x**3,
          lambda x: x**4]

    df =  [lambda x: np.ones(np.shape(x)),
           lambda x: 2*x,
           lambda x: 3*x**2,
           lambda x: 4*x**3]

    ddf =  [lambda x: np.zeros(np.shape(x)),
            lambda x: 2*np.ones(np.shape(x)),
            lambda x: 6*x,
            lambda x: 12*x**2]

    
    for i in xrange(len(f)):
        vs = UniformLagrangeVectorSpace(i+2)
        check_interpolation_derivative(vs, f[i], df[i])
        check_interpolation_derivative(vs, f[i], ddf[i], 2)

        vs *= 3
        check_interpolation_derivative(vs, f[i], df[i])
        check_interpolation_derivative(vs, f[i], ddf[i], 2)

        vs = UniformLagrangeVectorSpace(i+2)
        vs = AffineVectorSpace(vs, .3, .5)
        check_interpolation_derivative(vs, f[i], df[i])
        check_interpolation_derivative(vs, f[i], ddf[i], 2)
        
        
test_derivative()
