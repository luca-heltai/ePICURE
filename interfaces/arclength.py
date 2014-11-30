import sys
import numpy as np
import math
sys.path.append('../utilities')
from utilities.arc_lib import *


class ArcLengthParametrizer:
    """A python interface used to reparametrize a one dimensional curve. We
    assume that the curve parameter goes from 0 to 1 and we want to find the
    new curve with an arclength parametrization. The new parametrization will
    go again from 0 to 1 but with a constant velocity = L which is the overall
    length of the curve that we assume to be a constant.
    For the time-being we don't consider any time dependence. The curve is only
    space dependent. If you have something time dependent you are asked to query
    this class time by time.
    """
    def __init__(self, init_curve, arcfactor = 10):
		"""The constructor needs an initital curve and we ask for a factor that we will use to 
		numerically approximate the arclength, therefore we have called it arcfactor. We expect the 
		initial curve to have a field called knot_vector which is a numpy array."""
		self.init_curve = init_curve
		self.arcfactor = arcfactor
		self.nknot = self.init_curve.knot_vector.shape[0]
		self.new_knot_vector = np.zeros(self.nknot)

    def reparametrize(self):		
		"""This function compute the reparametrization. It queries the LS_assembler and LS_solver to compute
		the new knot vector that you can use to build an arclength reparametrized curve. """
		self.compute_arclength()
		self.reparametrization_LS_assembler()
		self.new_knot_vector = self.reparametrization_LS_solver()
		return self.new_knot_vector

    def compute_arclength(self):
		"""This function compute the overall length of the curve. We choose to do this
		in a very easy way. We consider a very large amount of points and we compute the
		length of the original curve considering it as piecewise rectilinear function.
		We return a numpy array holding the position in the old parametrization and the
		length of the curve at that point."""
		self.points_s = np.array([np.linspace(0, 1, self.arcfactor * self.arcfactor * self.nknot), 
								 np.zeros(self.arcfactor * self.arcfactor * self.nknot)])
		self.points_s = self.points_s.transpose()
		self.points_s[0,1] = 0.
		for i in range(1, self.points_s.shape[0]):
			P1 = self.init_curve.value(self.points_s[i - 1, 0])
			P2 = self.init_curve.value(self.points_s[i, 0])
			self.points_s[i,1] = np.linalg.norm(P2 - P1) + self.points_s[i-1,1]

		return self.points_s

    def find_s(self, sval):
		"""We call sval the arclength value we are quering. This function returns the tval, i.e.
		the value of the original parameter for which the desired arclength is achieved.
		We search over the array self.points_s for the closest arclength value to the desired 
		sval."""
		error = 1.
		for i in range(0, self.points_s.shape[0]):
			if(abs(self.points_s[i,1] - sval) < error):
				error = abs(self.points_s[i,1] - sval)
				tval = self.points_s[i,0]
		return tval

    def reparametrization_LS_assembler(self):
		"""In this function we compute the arclength reparametrization by mean of a Least Square 
		problem."""
		s_array = np.linspace(0, self.points_s[-1,1], self.arcfactor * self.nknot)
		self.point_ls = list()

		for i in range(0, s_array.shape[0]):
			tval = self.find_s(s_array[i])
			# The curve should have a value( or __call__ if prefered) method that we can query to know its value in space
			cval = self.init_curve.value(tval)
			# In point_ls we store the value in space, the s_array and the tval obtained
			self.point_ls.append([cval, s_array[i], tval])

		#self.point_ls = np.array(self.point_ls)
		# We compute the number of elements in the system rectangular matrix(Bmatrix), it will have dim*s_array.shape[0] rows and dim*nknot columns.
		# We want it to be rectangular because we are approximating its resilution so we search for something that solves the reparametrization in a 
		# least square sense.
		Bmatrixnumelem = self.init_curve.dim * s_array.shape[0] * self.nknot * self.init_curve.dim
		self.matrixB = np.zeros(Bmatrixnumelem).reshape(self.init_curve.dim * s_array.shape[0], self.nknot * self.init_curve.dim)
		self.rhsinit = np.zeros(self.init_curve.dim * s_array.shape[0])

		for i in range(0, s_array.shape[0]):
			for j in range(0, self.nknot):
				for d in range(0, self.init_curve.dim):
					self.matrixB[self.init_curve.dim * i + d, j + d * self.nknot] = self.init_curve.basis_functions[j].value(self.point_ls[i][2])

			for d in range(0, self.init_curve.dim):
				self.rhsinit[self.init_curve.dim * i + d] = self.point_ls[i][0][d]


    def reparametrization_LS_solver(self):
		"""In this method we solve the LS system assembled before. The result should be the new knot_vector that
		will form the arclength reparametrized curve."""
		res = np.linalg.lstsq(self.matrixB, self.rhsinit)
		new_knot_vector = res[0]
		return new_knot_vector

