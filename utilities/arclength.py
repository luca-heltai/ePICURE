import sys
import numpy as np
import math
from utilities.arc_lib import *
from utilities.matrices import *

class ArcLengthParametrizer(object):
    """A python utility used to reparametrize a one dimensional curve. We
    assume that the curve parameter goes from 0 to 1 and we want to find the
    new curve with an arclength parametrization. The new parametrization will
    go again from 0 to 1 but with a constant velocity = L which is the overall
    length of the curve that we assume to be a constant.
    For the time-being we don't consider any time dependence. The curve is only
    space dependent. If you have something time dependent you are asked to query
    this class time by time. We assume that the curve is parametrized from 0 to 1.
    """
    def __init__(self, vector_space, init_control_points, arcfactor = 10):
		"""The constructor needs an initital curve and we ask for a factor that we will use to 
		numerically approximate the arclength, therefore we have called it arcfactor. We expect the
		curve to be composed by a vector space and a control point array.
		"""
		self.init_control_points = np.asmatrix(init_control_points)
		self.vector_space = vector_space
		self.arcfactor = arcfactor
		self.n_dofs = self.vector_space.n_dofs
		self.dim = self.init_control_points[0].shape[1]
		self.new_control_points = np.zeros(self.n_dofs * self.dim)
		self.curve = self.vector_space.element(self.init_control_points)
		#self.curve_der = self.vector_space.element_der(self.init_control_points,1)

    def reparametrize(self):		
		"""This function compute the reparametrization. It queries the LS_assembler and LS_solver to compute
		the new knot vector that you can use to build an arclength reparametrized curve. """
		print "Starting the reparametrization"
		self.compute_arclength()
		print "Assembling the LS system"
		self.reparametrization_LS_assembler()
		print "Solving the system"
		self.new_control_points = self.reparametrization_LS_solver()
		print "Preparing the solution"
		
		new_cp = np.asmatrix(np.array([np.array([self.new_control_points[i*self.dim + j] for j in range(self.dim)]) for i in range(self.n_dofs)]))
		return new_cp

    def compute_arclength(self):
		"""This function compute the overall length of the curve. We choose to do this
		in a very easy way. We consider a very large amount of points and we compute the
		length of the original curve considering it as piecewise rectilinear function.
		We return a numpy array holding the position in the old parametrization and the
		length of the curve at that point."""
		self.points_s = np.matrix([np.linspace(0, 1, self.arcfactor * self.arcfactor * self.n_dofs), 
								 np.zeros(self.arcfactor * self.arcfactor * self.n_dofs)])
		self.points_s = self.points_s.transpose()
		self.points_s[0,1] = 0.
		all_points = self.curve(np.array(self.points_s[:,0]))
		all_points_diff = np.diff(all_points,1,0);
		#print all_points_diff.shape, all_points.shape
		for i in range(1, self.points_s.shape[0]):
			self.points_s[i,1] = np.linalg.norm(all_points_diff[i-1,:]) + self.points_s[i-1,1]
		return self.points_s

    def find_s(self, sval):
		"""We call sval the arclength value we are quering. This function returns the tval, i.e.
		the value of the original parameter for which the desired arclength is achieved.
		We search over the array self.points_s for the closest arclength value to the desired 
		sval."""
		error = 1.
		index = np.argmin(np.abs(self.points_s[:,1] - sval))
		return self.points_s[index,0]
		# for i in range(0, self.points_s.shape[0]):
		# 	if(abs(self.points_s[i,1] - sval) < error):
		# 		error = abs(self.points_s[i,1] - sval)
		# 		tval = self.points_s[i,0]
		# return tval

    def reparametrization_LS_assembler(self):
		"""In this function we compute the arclength reparametrization by mean of a Least Square 
		problem."""
		
		s_array = np.linspace(0, self.points_s[-1,1], self.arcfactor * self.n_dofs)
		#self.point_ls = list()
		self.point_ls = np.asmatrix(np.zeros([s_array.shape[0],self.dim + 2]))
		tval = np.zeros(s_array.shape[0])
		for i in range(0, s_array.shape[0]):
			tval[i] = self.find_s(s_array[i])
			#rint tval
			# The curve should have a value( or __call__ if prefered) method that we can query to know its value in space
		self.point_ls[:,0:self.dim] = self.curve(tval)
		#self.point_ls[:,self.dim] = s_array[:,np.newaxis]
		self.point_ls[:,self.dim] = tval[:,np.newaxis]
			# In point_ls we store the value in space, the s_array and the tval obtained
			#self.point_ls.append([cval, s_array[i], tval])

		#self.point_ls = np.array(self.point_ls)
		# We compute the number of elements in the system rectangular matrix(Bmatrix), it will have dim*s_array.shape[0] rows and dim*nknot columns.
		# We want it to be rectangular because we are approximating its resilution so we search for something that solves the reparametrization in a 
		# least square sense.
		Bmatrixnumelem = self.dim * s_array.shape[0] * self.n_dofs * self.dim
		self.matrixB = np.zeros(Bmatrixnumelem).reshape(self.dim * s_array.shape[0], self.n_dofs * self.dim)
		self.rhsinit = np.zeros(self.dim * s_array.shape[0])
		computation_matrix = interpolation_matrix(self.vector_space, tval)
		for i in range(0, s_array.shape[0]):
			for j in range(0, self.n_dofs):
				for d in range(0, self.dim):
					self.matrixB[self.dim * i + d, j * self.dim + d] = computation_matrix[i,j]#self.vector_space.basis(j)(self.point_ls[i, self.dim])

			for d in range(0, self.dim):
				self.rhsinit[self.dim * i + d] = self.point_ls[i,d]


    def reparametrization_LS_solver(self):
		"""In this method we solve the LS system assembled before. The result should be the new knot_vector that
		will form the arclength reparametrized curve."""
		res = np.linalg.lstsq(self.matrixB, self.rhsinit)
		return res[0]

