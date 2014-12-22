from interfaces import *
from utilities import *
import numpy as np

def test_derivative():
	n=2**3+1
	knots = np.linspace(0,1,n)
	# vs = UniformLagrangeVectorSpace(10)
	vs = IteratedVectorSpace(UniformLagrangeVectorSpace(3),knots)

	fun=lambda t: t**2
	s = np.linspace(0,1,vs.n_dofs)
	curve = fun(s)

	x    = vs.element(curve)
	dx   = vs.element_der(curve,1)

	times=np.linspace(0,1,10)
	eps=10e-11
	assert ((x(times) - times**2)<eps).all()
	assert (abs(dx(times) - 2*times)<eps).all()
