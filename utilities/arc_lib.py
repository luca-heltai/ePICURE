# In this library we try to create all the functions and objects needed to reparamitrized a curve using the arc-length parametrization

import numpy as np
import math
import scipy.special as sp 
from scipy.interpolate import lagrange 
from numpy.polynomial.chebyshev import chebgauss

# This class create and handles the Peano kernel and the Peano function
class PeanoHandler(object):
	def __init__(self,k):
		self.k = k
		#self.function = function



	def __call__(self, x_value, t_value):
		result = (x_value - t_value)**self.k #/ np.math.factorial(self.k)
		return result

	def anal_integral(self, a, b, t_value):
		result = ((b - t_value)**(self.k + 1) / (self.k + 1)) - ((a - t_value)**(self.k + 1) / (self.k + 1)) 
		return result

	def peano_kernel(self, a, b, t_value, quad_rule_w, quad_rule_q):
		p=PeanoHandler(self.k)
		numeric_integral = np.zeros(len(t_value))
		for j in range(len(t_value)):
			numeric_integral[j] = np.sum(np.array([quad_rule_w[i] * p(quad_rule_q[i], t_value[j]) for i in range(0, len(quad_rule_q))]))
		#print numeric_integral
		#numeric_integral = np.sum(np.array([self.value(quad_rule_q[i],t_value)*quad_rule_w[i] for i in range(0, len(quad_rule_q))]))
		result = self.anal_integral(a, b, t_value) - numeric_integral
		return result

	def peano_simpson(self, a, b, quad_rule, t_value):
		numeric_integral = np.sum([[quad_rule.quad_rule[i].w[j] * self.value(quad_rule.quad_rule[i].q[j], t_value) for j in range(0,len(quad_rule.quad_rule[i].w))] for i in range(0, len(quad_rule.quad_rule))])
		
		result = self.anal_integral(a, b, t_value) - numeric_integral
		return result

# In this class we deal with the creation of the Legendre polynomial using the Bonnet's formula
class LegendreBasis(object):
	def __init__(self,N):
		self.N = N + 1 
		self.coefficients = self.compute_coefficients()
		self.function = np.poly1d(self.coefficients)

	def __call__(self,x):
		return self.function(x)

	def compute_coefficients(self):
		if(self.N < 3):
			if(self.N == 1):
				coeff = [1];
			elif(self.N == 2):
				coeff = np.array([1, 0]);
			else:
				raise NotImplementedError
		else: 
			coeff1 = np.zeros(self.N)
			coeff2 = np.zeros(self.N)
			
			coeff1[0] = 1.
			coeff2[1] = 1.
			coeff = coeff2
			for i in range (3,self.N+1):
				coeff = np.zeros(self.N)
				#print coeff1, coeff2
				for j in range(1,i):
					#coeff[1:] += coeff2[j-1] * (2*i+1.)/(i+1) 
					coeff[j] += coeff2[j-1] * (2*i+1.)/(i+1) 
					#print coeff[j], coeff2[j-1], j-1
				for j in range(2,i):	
					coeff[j-2] += -coeff1[j-2] * (i)/(i+1)
				coeff1=coeff2
				coeff2=coeff
			coeff = coeff[::-1]

		return coeff

	def get_roots(self):
		return np.roots(self.coefficients)

# The lagrange basis handler
class lagrange_basis(object):
	def __init__(self,q):
		self.q=q
	def __call__(self,i):
		myf=self.q*0
		myf[i]=1
		return lagrange(self.q,myf)
	def integral(self,i,a,b):
		return np.polyint(self(i)) (b)-np.polyint(self(i)) (a)


class MyIntegral(object):
	def __init__(self,q,w):
		self.q = q
		self.w = w

	def __call__(self,f):
		return np.sum(self.w * f(self.q))

# The integral using the weights obtained with the lagrangian basis 
class MyLagIntegral(object):

	def __init__(self,q):
		self.q=q
		self.w=self.compute_weights()

	def compute_weights(self):
		li=lagrange_basis(self.q)
		return np.array([li.integral(i,0,1) for i in range(len(self.q))]) 


	def __call__(self,f):
		return np.sum(self.w*f(self.q)) # definition of approximated integral that we have done during lectures

# The same as bedore but rescaled from a to b
class MyLagScaledIntegral(MyLagIntegral):

	def __init__(self,q,a,b):
		#print q,b,a
		self.q=q
		self.w=self.compute_weights() * (b-a)
		if(b-a < 0):
			print "error!!!!!!!!!!!!", a, b
		self.q= q * (b-a) + a
		self.a=a
		self.b=b


	def get_weights(self):
		return self.w

	def get_quadrature_points(self):
		return self.q

# The iterated formula
class IterativeLagQuadrature(object):

	def __init__(self,a,b,Nq,N):
		self.a = a
		self.b = b
		self.intervals = np.linspace(a, b, N+1)
		self.Nq = Nq
		self.q = np.linspace(0,1,Nq)
		self.quad_rule = [MyLagScaledIntegral(self.q, self.intervals[i-1], self.intervals[i]) for i in range(1, len(self.intervals))]
		self.compute_weights()

	def __call__(self,f):
		return np.sum(self.ww * f(self.qq))
		#return np.sum([self.quad_rule[i].w * f(self.quad_rule[i].q) for i in range(0, len(self.quad_rule))])

	def compute_weights(self):
		self.ww=[]
		self.qq=[]
		for i in range(len(self.quad_rule)):
			self.ww.extend(self.quad_rule[i].w)
			self.qq.extend(self.quad_rule[i].q)
		self.ww = np.array(self.ww)
		self.qq = np.array(self.qq)

	def get_weights(self):
		return self.ww

	def get_quadrature_points(self):
		return self.qq

# The iterated formula using the Chebyshev poitns
class IterativeLagChebQuadrature(object):

	def __init__(self,a,b,Nq,N = 1):
		self.a = a
		self.b = b
		self.intervals = np.linspace(a,b,N+1)
		self.Nq = Nq
		self.q = chebgauss(N)[0] * 0.5 + 0.5
		self.quad_rule = [MyLagScaledIntegral(self.q, self.intervals[i-1], self.intervals[i]) for i in range(1, len(self.intervals))]
		self.compute_weights()

	def __call__(self,f):
		return np.sum(self.ww * f(self.qq))
		#return np.sum([self.quad_rule[i].w * f(self.quad_rule[i].q) for i in range(0, len(self.quad_rule))])

	def compute_weights(self):
		self.ww=[]
		self.qq=[]
		for i in range(len(self.quad_rule)):
			self.ww.extend(self.quad_rule[i].w)
			self.qq.extend(self.quad_rule[i].q)
		self.ww = np.array(self.ww)
		self.qq = np.array(self.qq)

	def get_weights(self):
		return self.ww

	def get_quadrature_points(self):
		return self.qq





class BasisFunction(object):
# Just the template for the BasisFunction, we consider a basic parametrization between 0 and 1, i.e. t in [0,1]. This function must have a value and a derivative.
	def __init__(self, knot_vector, my_num):
		self.my_num = np.array(my_num)
		self.knot_vector = np.array(knot_vector)


	def value(self, xval):
		return 1

	def derivative(self, xval):
		return 0

	def integral(self, xval):
		return 0

class BernsteinBasisFunction(BasisFunction):

	def value(self, xval):
		my_i = self.my_num
		ntot = self.knot_vector.shape[0] - 1
		return sp.binom(ntot, my_i) * xval**my_i * (1-xval)**(ntot - my_i)

	def gradient(self, xval, ntot, my_num):
		return sp.binom(ntot, my_i) *(my_i * xval**(my_i - 1) * 
			(1-xval)**(ntot - my_i) + (ntot - my_i) * xval**my_i * (1-xval)**(ntot - my_i -1))

class BernsteinBasisFunctionPowerLaw(BasisFunction):

	def value(self, xval):
		my_i = self.my_num
		ntot = self.knot_vector.shape[0] - 1
		bf = 0
		for j in range(my_i, ntot+1):
			bf += sp.binom(ntot, j) * sp.binom(j, my_i) * xval**(j) * (-1.)**(j-my_i) 
		return bf

	def gradient(self, xval, ntot, my_num):
		my_i = self.my_num
		ntot = self.knot_vector.shape[0] - 1
		dbf = 0
		for j in range(my_i, ntot+1):
			dbf += sp.binom(ntot, j) * sp.binom(j, my_i) * j * xval**(j-1) * (-1.)**(j-my_i) 
		return dbf


class BernsteinBasisFunctionRec(BasisFunction):

	def value(self, xval):
		my_i = self.my_num
		ntot = self.knot_vector.shape[0] - 1
		binm1 = sp.binom(ntot - 1, my_i) * xval**my_i * (1-xval)**(ntot - 1 - my_i)
		bim1nm1 = sp.binom(ntot - 1, my_i - 1) * xval**(my_i - 1) * (1-xval)**(ntot - 1 - (my_i - 1))
		return (1 - xval) * binm1 + xval * bim1nm1

	def gradient(self, xval, ntot, my_num):
		return sp.binom(ntot, my_i) *(my_i * xval**(my_i - 1) * 
			(1-xval)**(ntot - my_i) + (ntot - my_i) * xval**my_i * (1-xval)**(ntot - my_i -1))

	 
class LagrangianBasisFunction(BasisFunction):

	def value(self, xval):
		if self.my_num not in range(len(self.knot_vector)):
			raise NameError('Invali Lagrange Index')
		roots = self.knot_vector[range(self.my_num)+range(self.my_num + 1,len(self.knot_vector))]
		# we get the coefficients of the polynomial which has roots in roots
		wi = np.poly(roots)
		#we create the polynomial function
		pi = np.poly1d(wi)
		#finally we divide by its value in the point and we are done
		li = pi/pi(self.knot_vector[self.my_num])
		return li(xval)

	def derivative(self, xval):
		if self.my_num not in range(len(self.knot_vector)):
			raise NameError('Invali Lagrange Index')
		roots = self.knot_vector[range(self.my_num)+range(self.my_num + 1,len(self.knot_vector))]
		# we get the coefficients of the polynomial which has roots in roots
		wi = np.poly(roots)
		wid = np.polyder(wi)
		#we create the polynomial function
		pid = np.poly1d(wid)
		pi = np.poly1d(wi)
		#finally we divide by its value in the point and we are done
		lid = pid/pi(self.knot_vector[self.my_num])
		return lid(xval)





class MyCurve(object):
 	# the template for the Curve, we consider a basic parametrization between 0 and 1, i.e. t in [0,1]
	def __init__(self, basis_functions, control_points, dim):
		self.control_points = np.array(control_points)
		self.knot_vector = basis_functions[0].knot_vector
		self.basis_functions = basis_functions
		self.dim = dim

	def value(self, tval):
		val = 0 * self.control_points[0, :]
		for i in range(0, self.control_points.shape[0]):
			val = val + self.control_points[i,:] * self.basis_functions[i].value(tval)
		return val

	def derivative(self, tval):
		der = 0 * self.control_points[0, :]
		for i in range(0, self.control_points.shape[0]):
			der = der + self.control_points[i,:] * self.basis_functions[i].derivative(tval)
		return der




	





def Runge(xval):
	xval=np.array(xval)
	value =  1 / (1 + (10 * xval - 5)**2) 
	#value = [1 / (1 + (10 * xval[i] - 5)**2) for i in range(0, len(xval))]
	return value

def Linfnorm1d(x, yex, yapp):
	yex = np.array(yex)
	yapp = np.array(yapp)
	return np.max(abs(yex - yapp))



def ConvergenceError(error,h):
	assert(len(error)==len(h))
	#order=list()
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


		


