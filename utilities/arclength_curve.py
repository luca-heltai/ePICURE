import numpy as np
from numpy import linalg as LA
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from utilities.matrices import *
import interfaces
from utilities import least_square
from arclength import *
from scipy.interpolate import interp1d

class ArcLenghtCurve(object):
	"""A python utility used to construct and manipulate an arclength curve 
	in the space."""

	def __init__(self,vs):
		self.vs = vs
		self.space_dim = 0 #to be implemented
		self.gamma = ''
		self.coords = ''

	def from_lambda_to_coords(self, s_space=np.linspace(0,1,1025)):
		"""This method allows us to pass from lambda representation to
		coords representation of self.gamma."""
		assert (self.gamma != ''), \
				"self.gamma is not defined."
		# calculate the dimension of the space given the curve
		self.coords = least_square(self.gamma, self.vs, s_space)

	def from_coords_to_lambda(self):
		"""This method allows us to pass from coords representation to
		lambda representation of self.gamma."""
		assert (self.coords != ''), \
				"self.coords is not defined."
		self.gamma = self.vs.element(self.coords)

	def first_derivatives(self):
		""" If self.gamma is defined this method returns first, second, 
		and third derivative of the curve."""
		assert (self.gamma != ''), \
				"self.gamma is not defined."
		self.ds   = self.vs.element_der(self.coords,1)
		self.dds  = self.vs.element_der(self.coords,2)
		self.ddds = self.vs.element_der(self.coords,3)

	def first_derivative_modulus(self):
		assert (self.gamma != ''), \
			"self.gamma is not defined."
		self.first_derivatives()
		print ((self.ds(self.s_space).T).dot(self.ds(self.s_space))).diagonal()

	def first_normal_modulus(self):
		print ((self.normal(self.s_space)).dot(self.normal(self.s_space).T)).diagonal()

	def torsion(self):
		"""This method provides a lambda function representing the torsion
			of self.gamma."""
		self.first_derivatives()
		def tau(t):
			DX = self.ds(t)
			DDX = self.dds(t)
			DDDX = self.ddds(t)
			den = np.array(LA.norm(np.cross(DX.T,DDX.T).T,axis=0)**2)
			ids = np.array(den==0)
			den[ids]+=np.finfo(np.float64).eps
			return np.sum(np.cross(DX.T,DDX.T).T*DDDX,axis=0)/den
		self.tau = tau

	def curvature(self):
		"""This method provides a lambda function representing the curvature
			of self.gamma."""
		self.first_derivatives()
		def curv(t):
			DX = self.ds(t)
			DDX = self.dds(t)
			den = np.array((LA.norm(DX, axis=0))**3)
			ids = np.array(den==0)
			den[ids]+=np.finfo(np.float64).eps
			return LA.norm(np.cross(DX.T,DDX.T).T,axis=0)/den
		self.kappa = lambda t: curv(t)

	def reparametrise(self):
		"""This method reparametrise the curve in order to find an
		arclenght curve."""
		if(self.coords==''):
			self.from_lambda_to_coords()
		reparamCurve = ArcLengthParametrizer(self.vs,self.coords,\
							arcfactor = int(len(self.s_space))*100)
		self.gamma = lambda t: reparamCurve.curve(t)

	def curve_from_curvature(self,\
					start = 0, end = 1, \
					x0 = np.array([0,0,0]),\
					frenet0 = np.eye(3,3).reshape((-1,)),\
					L = 1):
		"""This method provides a curve starting from curvature and torsion."""
		M = lambda frenet,t: L*(np.array([\
		  		[0,         		self.kappa(t),   0			], \
				[-self.kappa(t), 	0,          	-self.tau(t)], \
				[0,         		self.tau(t),     0			]	] \
				).dot(frenet.reshape((3,3)))).reshape((-1,))
		self.frenet = odeint(M, frenet0, self.s_space).reshape((-1, 3,3))
		self.curve_discrete = x0 + \
			np.cumsum(self.frenet[:,:,0],axis=0)/(len(self.s_space))
		self.gamma = lambda t: np.array(interp1d(self.s_space,self.curve_discrete.T)(t)).T
		self.normal = lambda t: np.array(interp1d(self.s_space,self.frenet[:,:,1].T)(t)).T
		self.reparametrise()
		
	def plot(self, normal=False):
		"""This method plot slf.gamm in a 3D plot."""
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		plt.axis('equal')
		ax.plot(self.gamma(self.s_space)[0], self.gamma(self.s_space)[1], self.gamma(self.s_space)[2],'r', label='curve')
		if normal:
			ax.quiver(\
				self.gamma(self.s_space)[0]-self.normal(self.s_space)[:,0]*.01, \
				self.gamma(self.s_space)[1]-self.normal(self.s_space)[:,1]*.01, \
				self.gamma(self.s_space)[2]-self.normal(self.s_space)[:,2]*.01, \
				-self.normal(self.s_space)[:,0], -self.normal(self.s_space)[:,1], -self.normal(self.s_space)[:,2], \
				length=.01, label='normal', color='blue')
		ax.legend()
		plt.show()
		# plt.close()

	def LNorm(self):
		print np.sum(self.frenet[:,:,0]**2, axis=1)
		
	def NNorm(self):
		print np.sum(self.frenet[:,:,1]**2, axis=1)
		
	def BNorm(self):
		print np.sum(self.frenet[:,:,2]**2, axis=1)

class ALCFromCoords(ArcLenghtCurve):
	"""	Given a set of coordinates with respect to the basis
		of the vector space, this methot returns a lambda function
		representing the curve in the space."""
	def __init__(self, vs, coords):
		self.vs = vs
		self.coords = coords
		self.s_space = np.linspace(0,1,1025)
		self.from_coords_to_lambda()

class ALCFromLambda(ArcLenghtCurve):
	""" Given a vectorial lambda function this method returns the coordinates
	 	with respect to the basis of the vector space."""
	def __init__(self, vs, gamma, s_space=np.linspace(0,1,1025)):
		self.coords = ''
		self.vs = vs
		self.s_space=s_space
		self.gamma = lambda t : gamma(t).T
		self.from_lambda_to_coords()
		self.curvature()
		self.torsion()

class ALCFromKappaAndTau(ArcLenghtCurve):
	""" Given two lambda functions representing curvature and torsion,
	 	this class returns an arclenght curve. Moreover, curvature 
		and torsione are recalculated with respect the new curve."""
	def __init__(self, vs, kappa, tau,\
					x0 = np.array([0,0,0]),\
					frenet0 = np.eye(3,3).reshape((-1,)),\
					L = 1,\
					s_space = np.linspace(0,1,1025)):
		self.coords = ''
		self.kappa = kappa
		self.tau = tau
		self.vs = vs
		self.x0 = x0
		self.frenet0 = frenet0
		self.L = L
		self.s_space=s_space
		self.curve_from_curvature(
					start = 0, end = 1, \
					x0 = self.x0,\
					frenet0 = self.frenet0,\
					L = self.L)
		self.curvature()
		self.torsion()