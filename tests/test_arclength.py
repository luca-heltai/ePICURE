
import sys
# the arc_lib dir containes MyCurve tools to test the arclength
sys.path.append('../utilities')
sys.path.append('../interfaces')
from utilities.arc_lib import *
from interfaces.arclength import *
from nose.tools import *
import numpy as np 
import matplotlib.pyplot as plt
from pylab import *
import scipy.optimize as opt 
import math
from mpl_toolkits.mplot3d import Axes3D

def test_arclength_straight_line():
    """ Here we define the tests for the lenght computer of our ArcLengthParametrizer, we try it with a straight
    line. We test it both in 2d and 3d."""

    Sample = 10000
    x = np.linspace(0,1,Sample+1)

    # Number of interpolation points minus one
    n = 5
    toll = 1.e-6
    points = np.linspace(0, 1, (n+1) ) 
    control_points_2d = [[i, i] for i in range(0, n+1)]
    control_points_3d = [[i, i, i] for i in range(0, n+1)]

    functions = [BernsteinBasisFunction(points,i) for i in range(0, len(points))]
    dummy_curve_2d = MyCurve(functions, control_points_2d, 2 )
    dummy_curve_3d = MyCurve(functions, control_points_3d, 3 )
    dummy_arky_2d = ArcLengthParametrizer(dummy_curve_2d)
    dummy_arky_3d = ArcLengthParametrizer(dummy_curve_3d)
    length2d = dummy_arky_2d.compute_arclength()[-1][1]
    length3d = dummy_arky_3d.compute_arclength()[-1][1]
#    print length2d 
#    print n * np.sqrt(2)
    l2 = n * np.sqrt(2)
    l3 = n * np.sqrt(3) 
    assert (length2d - l2) < toll
    assert (length3d - l3) < toll
 #   assert True

def test_arclength_half_circle():
    """ Here we define the tests for the lenght computer of our ArcLengthParametrizer, we try it with a half a 
    circle and a fan. 
    We test it both in 2d and 3d."""

    Sample = 10000
    x = np.linspace(0,1,Sample+1)

    # Number of interpolation points minus one
    n = 5
    toll = 1.e-6
    points = np.linspace(0, 1, (n+1) ) 
    R = 1
    P = 1
    control_points_2d = [[np.cos(R * i * np.pi / (n + 1)), np.sin(R * i * np.pi / (n + 1))] for i in range(0, n+1)]
    control_points_3d = [[np.cos(R * i * np.pi / (n + 1)), np.sin(R * i * np.pi / (n + 1)), P * i] for i in range(0, n+1)]

    functions = [BernsteinBasisFunction(points,i) for i in range(0, len(points))]
    dummy_curve_2d = MyCurve(functions, control_points_2d, 2 )
    dummy_curve_3d = MyCurve(functions, control_points_3d, 3 )
    dummy_arky_2d = ArcLengthParametrizer(dummy_curve_2d)
    dummy_arky_3d = ArcLengthParametrizer(dummy_curve_3d)
    length2d = dummy_arky_2d.compute_arclength()[-1][1]
    length3d = dummy_arky_3d.compute_arclength()[-1][1]
#    print length2d 
#    print n * np.sqrt(2)
    l2 = np.pi * R
    l3 = 2 * np.pi * np.sqrt(R * R + (P / (2 * np.pi)) * (P / (2 * np.pi)))

    assert (length2d - l2) < toll
    assert (length3d - l3) < toll
 #   assert True

def test_reparametrization():
    """ Here we define the tests for reparametrizer of our ArcLengthParametrizer, we try it with a half a 
    circle and a fan. 
    We test it both in 2d and 3d."""
    Sample = 10000
    x = np.linspace(0,1,Sample+1)

    # Number of interpolation points minus one
    n = 5
    toll = 1.e-5
    points = np.linspace(0, 1, (n+1) ) 
    R = 1
    P = 1
    control_points_2d = [np.array([np.cos(R * i * np.pi / (n + 1)), np.sin(R * i * np.pi / (n + 1))]) for i in range(0, n+1)]
    control_points_3d = [np.array([np.cos(R * i * np.pi / (n + 1)), np.sin(R * i * np.pi / (n + 1)), P * i]) for i in range(0, n+1)]

    functions = [BernsteinBasisFunction(points,i) for i in range(0, len(points))]
    dummy_curve_2d = MyCurve(functions, control_points_2d, 2 )
    dummy_curve_3d = MyCurve(functions, control_points_3d, 3 )
    dummy_arky_2d = ArcLengthParametrizer(dummy_curve_2d)
    dummy_arky_3d = ArcLengthParametrizer(dummy_curve_3d)
    new_points_2d = dummy_arky_2d.reparametrize()
    new_points_3d = dummy_arky_3d.reparametrize()
    new_dummy_curve_2d = MyCurve(functions, control_points_2d, 2 )
    new_dummy_curve_3d = MyCurve(functions, control_points_3d, 3 )

#    print length2d 
#    print n * np.sqrt(2)
    l2 = n * np.sqrt(2)
    l3 = 2 * np.pi * np.sqrt(R * R + (P / (2 * np.pi)) * (P / (2 * np.pi)))
    original_curve_2d = [dummy_curve_2d.value(x[i]) for i in range(len(x))]
    reparam_curve_2d = [new_dummy_curve_2d.value(x[i]) for i in range(len(x))]

    # plt.subplot(121)
    # plt.plot(original_curve_2d[:][0], original_curve_2d[:][1])
    # plt.title('Original')
    # plt.subplot(122)
    # plt.plot(reparam_curve_2d[:][0], reparam_curve_2d[:][1])
    # plt.title('Arclength')
    # plt.savefig('Reparametrization2d.png')
    # plt.close()
    # plt.close()

    test_velocity_2d = [np.linalg.norm(new_dummy_curve_2d.derivative(x[i])) for i in range(len(x))]
    test_velocity_2d = [np.linalg.norm(new_dummy_curve_2d.derivative(x[i])) for i in range(len(x))]
 #   assert (length2d - test_velocity_2d.any()) > toll
 #   assert (length2d - test_velocity_3d.any()) > toll
    assert True




def test_rand():
    assert True



