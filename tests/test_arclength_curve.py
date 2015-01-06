from interfaces import *
from utilities import *
import numpy as np

def test_from_lambda_to_lambda():
	
	n=2**7+1
	knots = np.linspace(0,1,n)
	vector_space = IteratedVectorSpace(UniformLagrangeVectorSpace(5),knots)
	
	function = lambda t: np.array([np.sin(2*np.pi*t),np.cos(2*np.pi*t), t])

	curve1 = ALCFromLambda(vector_space,function)
	curve2 = ALCFromCoords(vector_space,curve1.coords)

	s_space=np.linspace(0,1,1025)
	assert (abs(function(s_space)-curve2.gamma(s_space))<10e-10).all(), \
	        'Expecting Failure! From lambda to lambda is not sharp enought'

def test_from_lambda_to_curvature():

	vs = UniformLagrangeVectorSpace(9)
	vector_space = AffineVectorSpace(vs)
	s_space = np.linspace(0,1,100)

	k=lambda t: 2*np.pi +0*t
	tau=lambda t: 0
	circ = ALCFromKappaAndTau(vector_space,k,tau)

	function=lambda t: circ.gamma(t)
	curve = ALCFromLambda(vector_space, function)
	newCurve = ALCFromKappaAndTau(vector_space,curve.kappa,curve.tau)

	assert (abs(newCurve.kappa(s_space)-k(s_space))<0.1).all(), \
        'Expecting Failure! Coputation of curvature form lambda is not sharp enought'
	assert (abs(newCurve.tau(s_space)-tau(s_space))<10e-10).all(), \
        'Expecting Failure! Coputation of torsion form lambda is not sharp enought'
	assert (abs(newCurve.gamma(s_space)-curve.gamma(s_space))<10e-2).all(), \
        'Expecting Failure! Coputation of the curve is not sharp enought'
		
def test_curvature():
	
	s_space=np.linspace(0,1,30)
	vs = UniformLagrangeVectorSpace(10)
	vector_space = AffineVectorSpace(vs)
	kappa = lambda t : t*0 + 1
	function = lambda t: np.array([np.sin(t),np.cos(t),0*t])
	curve = ALCFromLambda(vector_space, function)
	assert (abs(curve.kappa(s_space)-kappa(s_space))<10e-8).all(), \
	        'Expecting Failure! Computation of curvature is not sharp enought'
