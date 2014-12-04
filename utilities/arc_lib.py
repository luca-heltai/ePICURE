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
