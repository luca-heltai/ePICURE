# In this library we try to create all the functions and objects needed to reparamitrized a curve using the arc-length parametrization

import numpy as np
import math
import sys
#sys.path.append('../interfaces')
from interfaces.vector_space import *
from interfaces.lagrange_vector_space import *
from utilities.matrices import *

class Curve(object):
	"""A python utilities that handle a generic curve. As inputs it 
	requires a vector space and control points.
    """

 	
	def __init__(self, vector_space, control_points):
		"""We need to check that the input data are consistent"""
		assert(len(control_points) == vector_space.n_dofs)

		self.vector_space = AffineVectorSpace(vector_space, 0, 1)
		#self.vector_space = vector_space
		self.control_points = control_points
		self.n_dofs = vector_space.n_dofs


	def __call__(self):
		return self.vector_space.element(self.control_points)

	def curve_derivative(self):
		return self.vector_space.element_der(self.control_points) 

	def get_vector_space(self):
		return self.vector_space

	def get_vector_space(self):
		return self.vector_space

	def get_curve_n_dofs(self):
		return self.n_dofs

	def get_dimension(self):
		return len(self.control_points[0])






	





def Runge(xval):
	xval=np.array(xval)
	value =  1 / (1 + (10 * xval - 5)**2) 
	#value = [1 / (1 + (10 * xval[i] - 5)**2) for i in range(0, len(xval))]
	return value

def Linfnorm1d(x, yex, yapp):
	yex = np.array(yex)
	yapp = np.array(yapp)
	return np.max(abs(yex - yapp))

	error=np.array(error)
	h=np.array(h)
	order=[(np.log(error[i-1])/np.log(error[i])/ (np.log(h[i-1])/np.log(h[i]))) for i in range(1,len(error))]
	# for i in range(1,len(error)):
	# 	e1 = error[i-1]
	# 	e2 = error[i]
	# 	h1 = h[i-1]
	# 	h2 = h[i]
	# 	order.append(np.log(e1)/np.log(e2) / (np.ln(h1)/np.ln(h2)))
	order = np.array(order)
	return order


		


