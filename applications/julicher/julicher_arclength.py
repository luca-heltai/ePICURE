import numpy as np 
from interfaces import *
from utilities import *
import matplotlib.pyplot as plt
# We load the original control points exported from blender
cp_original = np.load('julicher_CPs.npy')
cp_original_2 = np.load('julicher_CPs_2.npy')
print cp_original.shape
if(cp_original.shape[0]==8):
	cp_original[1:-1,:,:]=cp_original[2:,:,:]
	#cp_original[-1,:,:]=0
	#cp_original=np.squeeze(cp_original)
	cp_original = np.delete(cp_original,len(cp_original)-1,0)
#print cp_original[-1,:,:]

#print cp_original_2

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

knots = np.array([0.000000, 0.000000, 0.000000, 0.000000, 0.307692, 0.384615, 0.461538, 0.538462, 0.615385, 0.692308, 1.000000, 1.000000, 1.000000, 1.000000]) 
knots_short = np.array([0.000000, 0.000000, 0.000000, 0.000000, 0.3 ,  0.45, 0.7, 1.000000, 1.000000, 1.000000, 1.000000])
print len(knots_short)
vsl = BsplineVectorSpace(p, knots)
vsl_short = BsplineVectorSpace(p, knots_short)
print vsl_short.n_dofs

arky = ArcLengthParametrizer(vsl, cp_original_2)
arky_short = ArcLengthParametrizer(vsl_short, cp_original)
cp_al = arky_short.reparametrize()

np.save("julicher_CPs_arclength",cp_al)
#np.save('julicher_CPs_arclength',cp_al)
#cp_al = np.load('julicher_CPs_arclength.npy')
# We may need to set some constraints
#cp_al[0,:,:]=cp_original[0,:,:]
print cp_al.shape
tt = np.linspace(0,1,128)
for i in np.array(range(14))*10:
	print i
	curve_orig = vsl_short.element(cp_original[:,i,:])
	curve_al = vsl_short.element(cp_al[:,i,:])
	plt.plot(cp_original[:,i,0], cp_original[:,i,1],'ro-')
	#plt.plot(cp_original_2[:,i,0], cp_original_2[:,i,1],'r*')
	plt.plot(curve_orig(tt)[0], curve_orig(tt)[1],'g-')
	plt.plot(cp_al[:,i,0],cp_al[:,i,1],'bo-')
	plt.plot(curve_al(tt)[0], curve_al(tt)[1],'k-')
	plt.savefig('jul_cp_'+str(i)+'.png')
	plt.close()
	#print np.amax(np.abs(cp_original[:,i,:]))
