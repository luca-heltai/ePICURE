import numpy as np 
from interfaces import *
from utilities import *
import matplotlib.pyplot as plt
# We load the original control points exported from blender
cp_original = np.load('julicher_CPs.npy')

#print cp_original

print cp_original.shape

#BSpline parameters
n = 5
p = 3
# Number of least square points
n_ls = 140

# Open knot vector
# knots = np.zeros(n+2*p)
# knots[p:-p] = np.linspace(0,1,n)
# knots[0:p] = 0
# knots[-p::] = 1
knots = np.array([0.000000, 0.000000, 0.000000, 0.000000, 0.307692, 0.5, 0.692308, 1.000000, 1.000000, 1.000000, 1.000000]) 
# 0.384615, 0.461538, 0.538462, 0.615385,
vsl = BsplineVectorSpace(p, knots)
print vsl.n_dofs

arky = ArcLengthParametrizer(vsl, cp_original)

#cp_al = arky.reparametrize()
#np.save('julicher_CPs_arclength',cp_al)
cp_al = np.load('julicher_CPs_arclength.npy')
# We may need to set some constraints
cp_al[0,:,:]=cp_original[0,:,:]
print cp_al.shape
tt = np.linspace(0,1,128)
for i in np.array(range(14))*10:
	print i
	curve_orig = vsl.element(cp_original[:,i,:])
	curve_al = vsl.element(cp_al[:,i,:])
	plt.plot(cp_original[:,i,0], cp_original[:,i,1],'ro-')
	plt.plot(curve_orig(tt)[0], curve_orig(tt)[1],'g-')
	plt.plot(cp_al[:,i,0],cp_al[:,i,1],'bo-')
	plt.plot(curve_al(tt)[0], curve_al(tt)[1],'k-')
	plt.savefig('jul_cp_'+str(i)+'.png')
	plt.close()
	print np.amax(np.abs(cp_original[:,i,:]))
