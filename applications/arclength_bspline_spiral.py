from interfaces import *
from utilities import *
from numpy import *
from numpy.linalg import lstsq
from matplotlib.pylab import *
from mpl_toolkits.mplot3d import Axes3D

# Spiral parameters:
nturns = 6.0
heigth = 1.0
radius = 1.0

# Spiral analytical expression
cx = lambda x: radius*exp(x)*sin(nturns*2*pi*x)
cy = lambda y: radius*2*exp(y)*cos(nturns*2*pi*y)
cz = lambda z: heigth*exp(z)*z

# BSpline parameters
n = 21
p = 3
# Number of least square points
n_ls = 140

# Open knot vector
knots = zeros(n+2*p)
knots[p:-p] = linspace(0,1,n)
knots[0:p] = 0
knots[-p::] = 1

vs = BsplineVectorSpace(p, knots)

# Least square parameter points
t = linspace(0,1,n_ls) 

# Least square points of the curve
F = array([cx(t), cy(t), cz(t)])

# Least square matrix
M = interpolation_matrix(vs, t)

# Control points and curve
CP = lstsq(M, F.T)[0]

arky = ArcLengthParametrizer(vs, CP)
CP_al = np.asarray(arky.reparametrize())
curve = vs.element(CP)
print (CP.shape)
curve_al = vs.element(CP_al)
arky_fixed = ArcLengthParametrizer(vs, CP, 10, 1)
CP_al_lf = np.asarray(arky_fixed.reparametrize())
new_arky_fixed = ArcLengthParametrizer(vs, CP_al_lf)
new_arky_fixed.reparametrize()
new_arky = ArcLengthParametrizer(vs, CP_al)
new_arky.reparametrize()
plt.plot(arky.points_s[:,0],arky.points_s[:,1],label='original')
plt.plot(new_arky.points_s[:,0],new_arky.points_s[:,1],label='reparametrized')
plt.plot(new_arky_fixed.points_s[:,0],new_arky_fixed.points_s[:,1],label='reparametrized_lf')
plt.legend()
plt.savefig('BSplinearclength.png')
plt.close()
plt.close()
#print (new_arky_fixed.points_s[-1,1], arky.points_s[-1,1])
print (CP.shape, type(CP))
print (CP_al.shape, type(CP_al))
# Approximated curve at points
C = curve(t)
C_al = curve(t)
# Plot curve, approximated curve and control polygon
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(C[0,:], C[1,:], C[2,:],'b')
ax.plot(F[0,:], F[1,:], F[2,:],'r')
ax.plot(CP[:,0], CP[:,1], CP[:,2], 'bo-')
#savefig('fig.png')
#plt.close()
#ax.plot(F[0,:], F[1,:], F[2,:])
ax.plot(CP_al[:,0], CP_al[:,1],CP_al[:,2], 'go-')
ax.plot(CP_al_lf[:,0], CP_al_lf[:,1],CP_al_lf[:,2], 'ko-')
savefig('fig_al.png')

