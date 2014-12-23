import sys
# the arc_lib dir containes MyCurve tools to test the arclength
sys.path.append('../utilities')
sys.path.append('../interfaces')
from utilities.arc_lib import *
from utilities.arclength import *
from nose.tools import *
import numpy as np 
import matplotlib.pyplot as plt
from pylab import *
import scipy.optimize as opt 
import math
from mpl_toolkits.mplot3d import Axes3D


# def test_arclength_half_circle():
#     """ Here we define the tests for the lenght computer of our ArcLengthParametrizer, we try it with a half a 
#     circle and a fan. 
#     We test it both in 2d and 3d."""


#     # Number of interpolation points minus one
#     n = 5
#     toll = 1.e-6
#     points = np.linspace(0, 1, (n+1) ) 
#     R = 1
#     P = 1
#     control_points_2d = np.asmatrix(np.zeros([n+1,2]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
#     control_points_2d[:,0] = np.transpose(np.matrix([R*np.cos(1 * i * np.pi / (n + 1))for i in range(n+1)]))
#     control_points_2d[:,1] = np.transpose(np.matrix([R*np.sin(1 * i * np.pi / (n + 1))for i in range(n+1)]))

#     control_points_3d = np.asmatrix(np.zeros([n+1,3]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
#     control_points_3d[:,0] = np.transpose(np.matrix([R*np.cos(1 * i * np.pi / (n + 1))for i in range(n+1)]))
#     control_points_3d[:,1] = np.transpose(np.matrix([R*np.sin(1 * i * np.pi / (n + 1))for i in range(n+1)]))
#     control_points_3d[:,2] = np.transpose(np.matrix([P*i for i in range(n+1)]))

#     vsl = AffineVectorSpace(UniformLagrangeVectorSpace(n+1),0,1)
#     dummy_arky_2d = ArcLengthParametrizer(vsl, control_points_2d)
#     dummy_arky_3d = ArcLengthParametrizer(vsl, control_points_3d)
#     length2d = dummy_arky_2d.compute_arclength()[-1,1]
#     length3d = dummy_arky_3d.compute_arclength()[-1,1]
# #    print length2d 
# #    print n * np.sqrt(2)
#     l2 = np.pi * R
#     l3 = 2 * np.pi * np.sqrt(R * R + (P / (2 * np.pi)) * (P / (2 * np.pi)))
#     print length2d, l2
#     print length3d, l3
#     assert (length2d - l2) < toll
#     assert (length3d - l3) < toll
# #   assert True

# def test_reparametrization():
#     """ Here we define the tests for reparametrizer of our ArcLengthParametrizer, we try it with a half a 
#     circle and a fan. 
#     We test it both in 2d and 3d."""
#     R = 1
#     P = 1
#     toll = 1.e-6

#     n = 10
#     ii = [0.5,0.8,4,6,7,8,9,9.6,10,10.1,11]
#     control_points_3d = np.asmatrix(np.zeros([n+1,3]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
#     control_points_3d[:,0] = np.transpose(np.matrix([R*np.cos(5*i * np.pi / (n + 1))for i in ii]))
#     control_points_3d[:,1] = np.transpose(np.matrix([R*np.sin(5*i * np.pi / (n + 1))for i in ii]))
#     control_points_3d[:,2] = np.transpose(np.matrix([P*i for i in range(n+1)]))
#     #control_points_3d[3,:] += 32
#     #print control_points_3d[0]
#     vsl = AffineVectorSpace(UniformLagrangeVectorSpace(n+1),0,1)
#     arky = ArcLengthParametrizer(vsl, control_points_3d)
#     new_control_points_3d = arky.reparametrize()

#     new_arky = ArcLengthParametrizer(vsl, new_control_points_3d)
#     new_new_control_points_3d = arky.reparametrize()
#     tt = np.linspace(0, 1, 128)

#     new_new_vals = vsl.element(new_new_control_points_3d)(tt)
#     #print vals
#     new_vals = vsl.element(new_control_points_3d)(tt)
#     #print vals.shape, new_vals.shape
#     assert np.amax(np.abs(new_new_vals-new_vals)) < toll

# #     Sample = 10000
# #     x = np.linspace(0,1,Sample+1)

# #     # Number of interpolation points minus one
# #     n = 5
# #     toll = 1.e-2
# #     points = np.linspace(0, 1, (n+1) ) 
# #     R = 1
# #     P = 1
# #     control_points_2d = [np.array([np.cos(R * i * np.pi / (n + 1)), np.sin(R * i * np.pi / (n + 1))]) for i in range(0, n+1)]
# #     control_points_3d = [np.array([np.cos(R * i * np.pi / (n + 1)), np.sin(R * i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]

# #     functions = [BernsteinBasisFunction(points,i) for i in range(0, len(points))]
# #     dummy_curve_2d = MyCurve(functions, control_points_2d, 2 )
# #     dummy_curve_3d = MyCurve(functions, control_points_3d, 3 )
# #     dummy_arky_2d = ArcLengthParametrizer(dummy_curve_2d)
# #     dummy_arky_3d = ArcLengthParametrizer(dummy_curve_3d)
# #     new_points_2d = dummy_arky_2d.reparametrize()
# #     new_points_3d = dummy_arky_3d.reparametrize()
# #     new_dummy_curve_2d = MyCurve(functions, control_points_2d, 2 )
# #     new_dummy_curve_3d = MyCurve(functions, control_points_3d, 3 )

# # #    print length2d 
# # #    print n * np.sqrt(2)
# #     l2 = n * np.sqrt(2)
# #     l3 = 2 * np.pi * np.sqrt(R * R + (P / (2 * np.pi)) * (P / (2 * np.pi)))
# #     original_curve_2d = [dummy_curve_2d.value(x[i]) for i in range(len(x))]
# #     reparam_curve_2d = [new_dummy_curve_2d.value(x[i]) for i in range(len(x))]

# #     # plt.subplot(121)
# #     # plt.plot(original_curve_2d[:][0], original_curve_2d[:][1])
# #     # plt.title('Original')
# #     # plt.subplot(122)
# #     # plt.plot(reparam_curve_2d[:][0], reparam_curve_2d[:][1])
# #     # plt.title('Arclength')
# #     # plt.savefig('Reparametrization2d.png')
# #     # plt.close()
# #     # plt.close()
# #     x=np.linspace(0.25,0.75,128)
# #     orig_velocity_2d = np.array([np.linalg.norm(dummy_curve_2d.derivative(x[i])) for i in range(len(x))])
# #     test_velocity_2d = np.array([np.linalg.norm(new_dummy_curve_2d.derivative(x[i])) for i in range(len(x))])
# #     test_velocity_3d = np.array([np.linalg.norm(new_dummy_curve_3d.derivative(x[i])) for i in range(len(x))])
# #     print orig_velocity_2d, test_velocity_2d

# #     #assert (l2 - test_velocity_2d.all()) < toll
# #     #assert (l3 - test_velocity_3d.all()) < toll
# #     assert True




# def test_known_parametrization():
#     R = 1
#     P = 1
#     toll = 1.e-6

#     n = 10
#     ii = np.linspace(0,1,n+1)
#     control_points_3d = np.asmatrix(np.zeros([n+1,3]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
#     control_points_3d[:,0] = np.transpose(np.matrix([R*np.cos(5*i * np.pi / (n + 1))for i in ii]))
#     control_points_3d[:,1] = np.transpose(np.matrix([R*np.sin(5*i * np.pi / (n + 1))for i in ii]))
#     control_points_3d[:,2] = np.transpose(np.matrix([P*i for i in range(n+1)]))
#     vsl = AffineVectorSpace(UniformLagrangeVectorSpace(n+1),0,1)
#     arky = ArcLengthParametrizer(vsl, control_points_3d)
#     new_control_points_3d = arky.reparametrize()

#     new_arky = ArcLengthParametrizer(vsl, new_control_points_3d)
#     new_new_control_points_3d = arky.reparametrize()
#     tt = np.linspace(0, 1, 128)

#     new_new_vals = vsl.element(new_new_control_points_3d)(tt)
#     #print vals
#     new_vals = vsl.element(new_control_points_3d)(tt)
#     #print vals.shape, new_vals.shape
#     assert np.amax(np.abs(new_new_vals-new_vals)) < toll

#     assert True



