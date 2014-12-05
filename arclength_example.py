import numpy as np
import math
import scipy.special as sp 
from scipy.interpolate import lagrange 
from numpy.polynomial.chebyshev import chebgauss
import sys
from utilities.arc_lib import *
from utilities.arclength import*
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



R = 1
P = 1

n = 18
ii = [(n+1) * np.cos(i * np.pi / (n + 1)) for i in range(n+1)]
control_points_3d = np.asmatrix(np.zeros([n+1,3]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
control_points_3d[:,0] = np.transpose(np.matrix([R*np.cos(5*i * np.pi / (n + 1))for i in ii]))
control_points_3d[:,1] = np.transpose(np.matrix([R*np.sin(5*i * np.pi / (n + 1))for i in ii]))
control_points_3d[:,2] = np.transpose(np.matrix([P*i for i in range(n+1)]))
#print control_points_3d[0]
vsl = IteratedVectorSpace(UniformLagrangeVectorSpace(3), np.linspace(0,1,10))
#vsl = AffineVectorSpace(UniformLagrangeVectorSpace(n+1),1,5)
arky = ArcLengthParametrizer(vsl, control_points_3d)
new_control_points_3d = arky.reparametrize()
new_arky = ArcLengthParametrizer(vsl, new_control_points_3d)
new_arky.reparametrize()
plt.plot(arky.points_s[:,0],arky.points_s[:,1],label='original')
plt.plot(new_arky.points_s[:,0],new_arky.points_s[:,1],label='reparametrized')
plt.legend()
plt.savefig('new_arclength.png')
plt.close()
plt.close()
print new_control_points_3d.shape
tt = np.linspace(0, 1, 128)
tt4 = 4 * tt + 1

vals = vsl.element(control_points_3d)(tt4)
#print vals
new_vals = vsl.element(new_control_points_3d)(tt4)
print vals.shape, new_vals.shape
x = np.squeeze(np.array(vals[:,0]))
y = np.squeeze(np.array(vals[:,1]))
z = np.squeeze(np.array(vals[:,2]))
new_x = np.squeeze(np.array(new_vals[:,0]))
new_y = np.squeeze(np.array(new_vals[:,1]))
new_z = np.squeeze(np.array(new_vals[:,2]))

print new_x.shape, x.shape, np.amax(np.abs(vals-new_vals),0)

#print control_points_3d[3]
#print x,y,z
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z,'r', label='test_curve')
ax.plot(np.squeeze(np.array(control_points_3d[:,0])),np.squeeze(np.array(control_points_3d[:,1])),np.squeeze(np.array(control_points_3d[:,2])),'r*',label='orig_cp')
ax.plot(new_x, new_y, new_z,'g', label='new_test_curve')
ax.plot(np.squeeze(np.array(new_control_points_3d[:,0])),np.squeeze(np.array(new_control_points_3d[:,1])),np.squeeze(np.array(new_control_points_3d[:,2])),'g*',label='new_cp')
ax.legend()
plt.savefig('test_curve.png')
plt.close()
plt.close()


