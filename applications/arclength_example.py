import numpy as np
import math
import scipy.special as sp 
from scipy.interpolate import lagrange 
from numpy.polynomial.chebyshev import chebgauss
import sys
from utilities import *
from interfaces import *
#from utilities.arclength import*
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



R = 1
P = 1
intervals=9
vs_order=2
n = (intervals*(vs_order)+1-1)
print (n)
ii = np.linspace(0,2,n+1)
#ii = [(n+1) * np.cos(i * np.pi / (n + 1)) for i in range(n+1)]
control_points_3d = np.asarray(np.zeros([n+1,2,3]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
control_points_3d[:,0,0] = np.array([R*np.cos(5*i * np.pi / (n + 1))for i in ii])
control_points_3d[:,0,1] = np.array([R*np.sin(5*i * np.pi / (n + 1))for i in ii])
control_points_3d[:,0,2] = np.array([P*i for i in range(n+1)])
control_points_3d[:,1,0] = np.array([R*np.cos(5*i * np.pi / (n + 1))for i in ii])
control_points_3d[:,1,1] = np.array([R*np.sin(5*i * np.pi / (n + 1))for i in ii])
control_points_3d[:,1,2] = np.array([2*P*i for i in range(n+1)])

#print control_points_3d[0]
vsl = IteratedVectorSpace(UniformLagrangeVectorSpace(vs_order+1), np.linspace(0,1,intervals+1))
print (vsl.n_dofs)
#vsl = AffineVectorSpace(UniformLagrangeVectorSpace(n+1),1,5)
#BSpline parameters
n = 17
p = 3
# Number of least square points
n_ls = 140

# Open knot vector
knots = np.zeros(n+2*p)
knots[p:-p] = np.linspace(0,1,n)
knots[0:p] = 0
knots[-p::] = 1

#vsl = BsplineVectorSpace(p, knots)
#print (vsl.n_dofs)
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
print (np.amax(np.abs(control_points_3d - new_control_points_3d )))
#print (np.squeeze(new_control_points_3d[:,0,:]))
tt = np.linspace(0, 1, 128)
tt4 = 4 * tt + 1
#print (tt4)
vals_1 = vsl.element(np.squeeze(control_points_3d[:,0,:]))(tt)
vals_2 = vsl.element(np.squeeze(control_points_3d[:,1,:]))(tt)
new_vals_1 = vsl.element(np.squeeze(new_control_points_3d[:,0,:]))(tt)
new_vals_2 = vsl.element(np.squeeze(new_control_points_3d[:,1,:]))(tt)

#print (vsl.element(np.squeeze(control_points_3d[:,1,:]))(tt4) == arky.curve(tt4)), #vals_1

#print (vals.shape, new_vals.shape)
x = np.squeeze(np.array(vals_1[0,:]))
y = np.squeeze(np.array(vals_1[1,:]))
z = np.squeeze(np.array(vals_1[2,:]))
#print (x.shape)
new_x = np.squeeze(np.array(new_vals_1[0,:]))
new_y = np.squeeze(np.array(new_vals_1[1,:]))
new_z = np.squeeze(np.array(new_vals_1[2,:]))

#print (new_x.shape, x.shape, np.amax(np.abs(vals-new_vals),0))

#print (control_points_3d[3])
#print (x,y,z)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z,'r', label='test_curve')
ax.plot(np.squeeze(np.array(control_points_3d[:,0,0])),np.squeeze(np.array(control_points_3d[:,0,1])),np.squeeze(np.array(control_points_3d[:,0,2])),'r*-',label='orig_cp')
ax.plot(new_x, new_y, new_z,'g', label='new_test_curve')
ax.plot(np.squeeze(np.array(new_control_points_3d[:,0,0])),np.squeeze(np.array(new_control_points_3d[:,0,1])),np.squeeze(np.array(new_control_points_3d[:,0,2])),'g*-',label='new_cp')
ax.legend()
plt.savefig('test_curve_1.png')
plt.close()
plt.close()

#print (vals_1.shape, new_vals_1.shape)
x = np.squeeze(np.array(vals_2[0,:]))
y = np.squeeze(np.array(vals_2[1,:]))
z = np.squeeze(np.array(vals_2[2,:]))
new_x = np.squeeze(np.array(new_vals_2[0,:]))
new_y = np.squeeze(np.array(new_vals_2[1,:]))
new_z = np.squeeze(np.array(new_vals_2[2,:]))
fig = plt.figure()
#print (new_x)
ax = fig.gca(projection='3d')
ax.plot(x, y, z,'r', label='test_curve')
ax.plot(np.squeeze(np.array(control_points_3d[:,1,0])),np.squeeze(np.array(control_points_3d[:,1,1])),np.squeeze(np.array(control_points_3d[:,1,2])),'r*-',label='orig_cp')
ax.plot(new_x, new_y, new_z,'g', label='new_test_curve')
ax.plot(np.squeeze(np.array(new_control_points_3d[:,1,0])),np.squeeze(np.array(new_control_points_3d[:,1,1])),np.squeeze(np.array(new_control_points_3d[:,1,2])),'g*-',label='new_cp')
ax.legend()
plt.savefig('test_curve_2.png')
plt.close()
plt.close()



