from __future__ import print_function
from numpy import ones, zeros, array, asarray
from scipy.misc import comb, derivative

# no xrange on python3
try:
	x_range = xrange
except NameError:
	x_range = range

#======================================================================#
# **Bernstein Basis**

def bernstein_poly(n,i,t):
	"""
	returns the `ith` Bernstein polynomial of degree n as a function of t
		B_(i,n)(t) = (^n _i)*t^i*(1-t)^(n-i)
	for i = 0, 1, ..., n, where
		(^n _i) = n! / (i!(n-i)!)
	There are n+1 nth-degree Bernstein polynomials.
	if i < 0 or i > n, B_(i,n) = 0
	"""
	if i < 0 or i > n:
		return zeros(len(t),)
	if n == 0:
		return ones(len(t),)
	return comb(n,i) * (t**i) * (1-t)**(n-i)


def bernstein_get_poly(n,k):
	"""
	returns a Bernstein polynomial (as a callable object)
	usage example:
		n = 3
		x = np.linspace(0,1,1025)
		y = [bernstein_get_poly(n,i)(x) for i in range(n+1)]
		plt.plot(x,y)
	"""
	def _bpoly(x):
		return bernstein_poly(n,k,x)
	return _bpoly


def bernstein_poly_r(n,k,t):
	"""
	Recursive Definition of the Bernstein Polynomials.
	The Bernstein polynomials of degree `n` can be defined by blending
	together two Bernstein polynomials of degree `n-1`.
	Returns the `kth` nth-degree Bernstein polynomial as:
		B_(k,n)(t) = (1-t)*B_(k,n-1)(t) + t*B(k-1,n-1)(t)
	"""
	if k < 0 or k > n:
		return zeros(len(t),)
	if n == 0 and k == 0:
		return ones(len(t),)
	return (1-t) * bernstein_poly_r(n-1,k,t) + t * bernstein_poly_r(n-1,k-1,t)


def bernstein_poly_pbasis(n,k,t):
	"""
	returns the `kth` Bernstein polynomial in terms of the power basis:
		B_(k,n)(t) = \sum_{i=k}^n (-1)^(i-k)*(^n _i)*(^i _k)*t^i
	"""
	if k < 0 or k > n:
		return zeros(len(t),)
	if n == 0:
		return ones(len(t),)
	return sum([((-1)**(i-k)) * comb(n,i) * comb(i,k) * (t**i) for i in x_range(k,n+1)])


def bernstein_derivative(n,i,x,h):
	"""
	returns the derivative of the i-th Bernstein basis, evaluated at `x`
	"""
	f = bernstein_get_poly(n,i)
	return derivative(f, x, dx=h, n=1, args=(), order=3)


def bernstein_approx(n,f,x):
	"""
	returns the Bernstein approximation of a continous function as:
		B_n(f, x) := \sum_{i=0}^n B(i,n,x) * f(i/n)
	"""
	return sum([bernstein_poly(n,i,x)*f(float(i)/n) for i in x_range(n)])



if __name__ == "__main__":
	# test case if invoked directly

	from numpy import linspace

	points = 1025
	n = 3
	x = linspace(0,1,points)

	print(repr(x))

	bp   = [bernstein_poly(n,i,x) for i in x_range(n+1)]
	bp_b = [bernstein_get_poly(n,i)(x) for i in x_range(n+1)]
	bp_p = [bernstein_poly_pbasis(n,i,x) for i in x_range(n+1)]
	bp_r = [bernstein_poly_r(n,i,x) for i in x_range(n+1)]

	print(repr(bp))
	print(repr(bp_b))
	print(repr(bp_p))
	print(repr(bp_r))

	n = 5
	i = 3
	x = 1
	h = 1e-6
	d = bernstein_derivative(n,i,x,h)

	print(repr(d))

	#bernstein_approx(n,f,x)

#EOF
