import sys
from interfaces import *
from utilities import *
from nose.tools import *
import numpy as np 
import matplotlib.pyplot as plt
from pylab import *
import scipy.optimize as opt 
import math
from mpl_toolkits.mplot3d import Axes3D


def test_arclength_half_circle():
    """ Here we define the tests for the lenght computer of our ArcLengthParametrizer, we try it with a half a 
    circle and a fan. 
    We test it both in 2d and 3d."""


    # Number of interpolation points minus one
    n = 5
    toll = 1.e-6
    points = np.linspace(0, 1, (n+1) ) 
    R = 1
    P = 1
    control_points_2d = np.asmatrix(np.zeros([n+1,2]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
    control_points_2d[:,0] = np.transpose(np.matrix([R*np.cos(1 * i * np.pi / (n + 1))for i in range(n+1)]))
    control_points_2d[:,1] = np.transpose(np.matrix([R*np.sin(1 * i * np.pi / (n + 1))for i in range(n+1)]))

    control_points_3d = np.asmatrix(np.zeros([n+1,3]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
    control_points_3d[:,0] = np.transpose(np.matrix([R*np.cos(1 * i * np.pi / (n + 1))for i in range(n+1)]))
    control_points_3d[:,1] = np.transpose(np.matrix([R*np.sin(1 * i * np.pi / (n + 1))for i in range(n+1)]))
    control_points_3d[:,2] = np.transpose(np.matrix([P*i for i in range(n+1)]))

    vsl = AffineVectorSpace(UniformLagrangeVectorSpace(n+1),0,1)
    dummy_arky_2d = ArcLengthParametrizer(vsl, control_points_2d)
    dummy_arky_3d = ArcLengthParametrizer(vsl, control_points_3d)
    length2d = dummy_arky_2d.compute_arclength()[-1,1]
    length3d = dummy_arky_3d.compute_arclength()[-1,1]
#    print (length2d)
#    print (n * np.sqrt(2))
    l2 = np.pi * R
    l3 = 2 * np.pi * np.sqrt(R * R + (P / (2 * np.pi)) * (P / (2 * np.pi)))
    print (length2d, l2)
    print (length3d, l3)
    assert (length2d - l2) < toll
    assert (length3d - l3) < toll
#   assert True

def test_reparametrization():
    """ Here we define the tests for reparametrizer of our ArcLengthParametrizer, we try it with a half a 
    circle and a fan. 
    We test it both in 2d and 3d."""
    R = 1
    P = 1
    toll = 1.e-6

    n = 10
    ii = [0.5,0.8,4,6,7,8,9,9.6,10,10.1,11]
    control_points_3d = np.asmatrix(np.zeros([n+1,3]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
    control_points_3d[:,0] = np.transpose(np.matrix([R*np.cos(5*i * np.pi / (n + 1))for i in ii]))
    control_points_3d[:,1] = np.transpose(np.matrix([R*np.sin(5*i * np.pi / (n + 1))for i in ii]))
    control_points_3d[:,2] = np.transpose(np.matrix([P*i for i in range(n+1)]))
    #control_points_3d[3,:] += 32
    #print control_points_3d[0]
    vsl = AffineVectorSpace(UniformLagrangeVectorSpace(n+1),0,1)
    arky = ArcLengthParametrizer(vsl, control_points_3d)
    new_control_points_3d = arky.reparametrize()

    new_arky = ArcLengthParametrizer(vsl, new_control_points_3d)
    new_new_control_points_3d = arky.reparametrize()
    tt = np.linspace(0, 1, 128)

    new_new_vals = vsl.element(new_new_control_points_3d)(tt)
    #print vals
    new_vals = vsl.element(new_control_points_3d)(tt)
    #print vals.shape, new_vals.shape
    assert np.amax(np.abs(new_new_vals-new_vals)) < toll

#     Sample = 10000
#     x = np.linspace(0,1,Sample+1)

#     # Number of interpolation points minus one
#     n = 5
#     toll = 1.e-2
#     points = np.linspace(0, 1, (n+1) ) 
#     R = 1
#     P = 1
#     control_points_2d = [np.array([np.cos(R * i * np.pi / (n + 1)), np.sin(R * i * np.pi / (n + 1))]) for i in range(0, n+1)]
#     control_points_3d = [np.array([np.cos(R * i * np.pi / (n + 1)), np.sin(R * i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]

#     functions = [BernsteinBasisFunction(points,i) for i in range(0, len(points))]
#     dummy_curve_2d = MyCurve(functions, control_points_2d, 2 )
#     dummy_curve_3d = MyCurve(functions, control_points_3d, 3 )
#     dummy_arky_2d = ArcLengthParametrizer(dummy_curve_2d)
#     dummy_arky_3d = ArcLengthParametrizer(dummy_curve_3d)
#     new_points_2d = dummy_arky_2d.reparametrize()
#     new_points_3d = dummy_arky_3d.reparametrize()
#     new_dummy_curve_2d = MyCurve(functions, control_points_2d, 2 )
#     new_dummy_curve_3d = MyCurve(functions, control_points_3d, 3 )

# #    print length2d 
# #    print n * np.sqrt(2)
#     l2 = n * np.sqrt(2)
#     l3 = 2 * np.pi * np.sqrt(R * R + (P / (2 * np.pi)) * (P / (2 * np.pi)))
#     original_curve_2d = [dummy_curve_2d.value(x[i]) for i in range(len(x))]
#     reparam_curve_2d = [new_dummy_curve_2d.value(x[i]) for i in range(len(x))]

#     # plt.subplot(121)
#     # plt.plot(original_curve_2d[:][0], original_curve_2d[:][1])
#     # plt.title('Original')
#     # plt.subplot(122)
#     # plt.plot(reparam_curve_2d[:][0], reparam_curve_2d[:][1])
#     # plt.title('Arclength')
#     # plt.savefig('Reparametrization2d.png')
#     # plt.close()
#     # plt.close()
#     x=np.linspace(0.25,0.75,128)
#     orig_velocity_2d = np.array([np.linalg.norm(dummy_curve_2d.derivative(x[i])) for i in range(len(x))])
#     test_velocity_2d = np.array([np.linalg.norm(new_dummy_curve_2d.derivative(x[i])) for i in range(len(x))])
#     test_velocity_3d = np.array([np.linalg.norm(new_dummy_curve_3d.derivative(x[i])) for i in range(len(x))])
#     print orig_velocity_2d, test_velocity_2d

#     #assert (l2 - test_velocity_2d.all()) < toll
#     #assert (l3 - test_velocity_3d.all()) < toll
#     assert True




def test_known_parametrization():
    R = 1
    P = 1
    toll = 2.e-3

    n = 10
    ii = np.linspace(0,1,n+1)
    control_points_3d = np.asarray(np.zeros([n+1,3]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
    print (control_points_3d.shape)
    control_points_3d[:,0] = np.array([R*np.cos(5*i * np.pi / (n + 1))for i in ii])
    control_points_3d[:,1] = np.array([R*np.sin(5*i * np.pi / (n + 1))for i in ii])
    control_points_3d[:,2] = np.array([P*i for i in range(n+1)])
    vsl = AffineVectorSpace(UniformLagrangeVectorSpace(n+1),0,1)
    arky = ArcLengthParametrizer(vsl, control_points_3d)
    new_control_points_3d = arky.reparametrize()

    #new_arky = ArcLengthParametrizer(vsl, new_control_points_3d)
    #new_new_control_points_3d = arky.reparametrize()
    tt = np.linspace(0, 1, 128)

    vals = vsl.element(control_points_3d)(tt)
    #print vals
    new_vals = vsl.element(new_control_points_3d)(tt)
    #print vals.shape, new_vals.shape
    print (np.amax((np.abs(vals-new_vals))))
    assert np.amax(np.abs(control_points_3d-new_control_points_3d))/P < toll

    #assert True

def test_more_known_parametrization_together():
    R = 1
    P = 1
    toll = 7.e-3
    intervals = 5
    vs_order = 2
    n = (intervals*(vs_order)+1-1)

    #n = 18
    ii = np.linspace(0,1,n+1)
    n_1 = 2
    n_2 = 4
    control_points_3d = np.asarray(np.zeros([n+1,n_1,n_2,3]))#[np.array([R*np.cos(5*i * np.pi / (n + 1)), R*np.sin(5*i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]
    for k in range(n_1):
        for j in range(n_2):
            control_points_3d[:,k,j,0] = np.array([R*np.cos(5*i * np.pi / (n + 1))for i in ii])
            control_points_3d[:,k,j,1] = np.array([R*np.sin(5*i * np.pi / (n + 1))for i in ii])
            control_points_3d[:,k,j,2] = np.array([(k+j+1)*P*i for i in range(n+1)])
    #vsl = IteratedVectorSpace(UniformLagrangeVectorSpace(vs_order+1), np.linspace(0,1,intervals+1))
    vsl = AffineVectorSpace(UniformLagrangeVectorSpace(n+1),0,1)
    arky = ArcLengthParametrizer(vsl, control_points_3d)
    new_control_points_3d = arky.reparametrize()

    #print control_points_3d.shape, new_control_points_3d.shape
    tt = np.linspace(0,1,128)
    for k in range(n_1):
        for j in range(n_2):
            vals = vsl.element(control_points_3d)(tt)
            new_vals = vsl.element(new_control_points_3d)(tt)
            print (np.amax(np.abs(vals-new_vals))/(k+j+1)/P, (k+j+1))
            assert np.amax(np.abs(vals-new_vals))/(k+j+1)/P < toll

def test_length_constraint():
    toll = 1e-5
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
    M = InterpolationMatrix(vs, t)

    # Control points and curve
    CP = lstsq(M, F.T, rcond=-1)[0]
    arky_fixed = ArcLengthParametrizer(vs, CP, 1)
    CP_al_lf = np.asarray(arky_fixed.reparametrize())
    new_arky_fixed = ArcLengthParametrizer(vs, CP_al_lf)
    new_arky_fixed.reparametrize()
    assert(np.abs(arky_fixed.lengths[0]-new_arky_fixed.lengths[0]) < toll)


def test_more_length_constraints():
    toll = 1e-5
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
    M = InterpolationMatrix(vs, t)

    # Control points and curve
    CP = lstsq(M, F.T, rcond=None)[0]
    CP2 = np.empty((CP.shape[0],2,CP.shape[1]))
    CP2[:,0,:] = CP
    CP2[:,1,:] = 2*CP
    arky_fixed1 = ArcLengthParametrizer(vs, CP2, 1)
    arky_fixed2 = ArcLengthParametrizer(vs, CP2, 2)
    CP_al_lf1 = np.asarray(arky_fixed1.reparametrize())
    CP_al_lf2 = np.asarray(arky_fixed2.reparametrize())
    new_arky_fixed1 = ArcLengthParametrizer(vs, CP_al_lf1)
    new_arky_fixed1.reparametrize()
    new_arky_fixed2 = ArcLengthParametrizer(vs, CP_al_lf2)
    new_arky_fixed2.reparametrize()
    assert(np.abs(arky_fixed1.lengths[0]-new_arky_fixed1.lengths[0]) < toll)
    assert(np.abs(arky_fixed1.lengths[1]-new_arky_fixed1.lengths[1]) < toll)
    print (arky_fixed1.lengths[1], new_arky_fixed2.lengths[0])
    assert(np.abs(arky_fixed1.lengths[1]-new_arky_fixed2.lengths[0]) < toll)
    assert(np.abs(arky_fixed1.lengths[1]-new_arky_fixed2.lengths[1]) < toll)



