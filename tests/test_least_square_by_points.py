from interfaces import *
from utilities import *
from numpy import *
from numpy.linalg import lstsq
from matplotlib.pylab import *
from mpl_toolkits.mplot3d import Axes3D

def test_ls_points():
    # Spiral parameters:
    toll = 1e-6
    nturns = 6.0
    heigth = 1.0
    radius = 1.0

    # Spiral analytical expression
    cx = lambda x: radius*sin(nturns*2*pi*x)
    cy = lambda y: radius*cos(nturns*2*pi*y)
    cz = lambda z: heigth*z

    # BSpline parameters
    n = 21
    p = 3
    # Number of least square points
    n_ls = 140

    # Open knot vector
    knots = r_[p*[0], linspace(0,1,n), p*[1]]

    vs = BsplineVectorSpace(p, knots)

    # Least square parameter points
    t = linspace(0,1,n_ls)

    # Least square points of the curve
    F = array([cx(t), cy(t), cz(t)])

    # Least square matrix
    M = interpolation_matrix(vs, t)

    # Control points and curve
    CP = lstsq(M, F.T, rcond=-1)[0]
    CP2 = least_square_by_points(vs, F.T, t)
    assert(np.linalg.norm(CP-CP2)<toll)