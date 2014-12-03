import math
from _Bas import bspkntins, bspdegelev, bspbezdecom, bspeval # Lowlevel Nurbs functions
from Util import  scale, translate, rotz, NURBSError

dependencies = '''This module requires:
	Numeric Python (NumPy)
'''

try:
    import numpy as np
except ImportError, value:
	print dependencies
	raise

class Crv(object):
    '''Construct a NURB curve and check the format.
    
 The NURB curve is represented by a 4 dimensional b-spline.

 INPUT:

    cntrl  - Control points, homogeneous coordinates (wx,wy,wz,w)
            [dim,nu] matrix
            dim is the dimension valid options are:
            2 .... (x,y)        2D cartesian coordinates
            3 .... (x,y,z)      3D cartesian coordinates   
            4 .... (wx,wy,wz,w) 4D homogeneous coordinates

    uknots - Knot sequence along the parametric u direction.

 NOTES:

    Its assumed that the input knot sequences span the
    interval [0.0,1.0] and are clamped to the control
    points at the end by a knot multiplicity equal to
    the spline order.'''

    def __init__(self, cntrl, uknots):

        # Placeholder for Bezier representation of self
        self._bezier = None

        # Force the u knot sequence to be a vector in ascending order
        # and normalise between [0.0,1.0]
        uknots = np.sort(np.asarray(uknots, np.float))
        nku = uknots.shape[0]
        uknots = (uknots - uknots[0])/(uknots[-1] - uknots[0])

        # TODO: insert more knot checks
        if uknots[0] == uknots[-1]:
            raise NURBSError, 'Illegal uknots sequence'

        self.uknots = uknots

        # Complain about control points that are less than 2D or more
        # than 4D. Fill 3rd dimension with zeros and 4th dimension
        # with ones if necessary.
        cntrl = np.asarray(cntrl, np.float)
        (dim, nu) = cntrl.shape
        if dim < 2 or dim > 4:
            raise NURBSError, 'Illegal control point format'
        elif dim < 4:
            self.cntrl = np.zeros((4, nu), np.float)
            self.cntrl[0:dim,:] = cntrl
            self.cntrl[-1,:] = np.ones((nu,))
        else:
            self.cntrl = cntrl

        # Spline degree
        self.degree = nku - nu - 1
        if self.degree < 0:
            raise NURBSError, 'NURBS order must be a positive integer'

    def trans(self, mat):
        "Apply the 4D transform matrix to the NURB control points."
        self.cntrl = np.dot(mat, self.cntrl)

    def reverse(self):
        "Reverse evaluation direction"
        self.cntrl = self.cntrl[:,::-1]
        self.uknots = 1 - self.uknots[::-1]
        
    def kntins(self, uknots):
        """Insert new knots into the curve
	NOTE: No knot multiplicity will be increased beyond the order of the spline"""
        if len(uknots):
            uknots = np.sort(np.asarray(uknots, np.float))
            if np.any(uknots < 0.) or np.any(uknots > 1.):
                raise NURBSError, 'NURBS curve parameter out of range [0,1]'
            self.cntrl, self.uknots = bspkntins(self.degree, self.cntrl, self.uknots, uknots)

    def degelev(self, degree):
        "Degree elevate the curve"
        if degree < 0:
            raise NURBSError, 'degree must be a positive number'
        if degree > 0:
            cntrl, uknots, nh = bspdegelev(self.degree, self.cntrl, self.uknots, degree)
            self.cntrl = cntrl[:,:nh + 1]
            self.uknots = uknots[:nh + self.degree + degree + 2]
            self.degree += degree

    def bezier(self, update = None):
        "Decompose curve to bezier segments and return overlaping control points"
        if update or not self._bezier:
            self._bezier = bspbezdecom(self.degree, self.cntrl, self.uknots)
        return self._bezier

    def bounds(self):
        "Return the boundingbox for the curve"
        ww = np.resize(self.cntrl[-1,:], (3, self.cntrl.shape[1]))
        cntrl = np.sort(self.cntrl[0:3,:]/ww)
        return np.asarray([cntrl[0,0], cntrl[1,0], cntrl[2,0],
                                cntrl[0,-1], cntrl[1,-1], cntrl[2,-1]], np.float)
                                
    def pnt3D(self, ut):
        "Evaluate parametric point[s] and return 3D cartesian coordinate[s]"
        val = self.pnt4D(ut)
        return val[:3]/val[3]

    def pnt4D(self, ut):
        "Evaluate parametric point[s] and return 4D homogeneous coordinates"
        ut = np.asarray(ut, np.float)
		
        if np.any(ut < 0.) or np.any(ut > 1.):
            raise NURBSError, 'NURBS curve parameter out of range [0,1]'
        return bspeval(self.degree, self.cntrl, self.uknots, ut)
                
    def plot(self, n = 25):
        """A simple plotting function for debugging purpose
	n = number of subdivisions.
	Depends on the matplotlib plotting library."""

        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plt

        pnts = self.pnt3D(np.arange(n + 1, dtype = np.float64)/n)
        knot = self.pnt3D(self.uknots)
        ctrl = self.cntrl[:3]/self.cntrl[3]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.set_title( "b-spline Curve, degree={0}".format(self.degree) )
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        ax.plot(pnts[0], pnts[1], pnts[2], label='parametric bspline')
        ax.plot(ctrl[0], ctrl[1], ctrl[2], 'ro-', label='control pts', linewidth=.5)
        ax.plot(knot[0], knot[1], knot[2], 'y+', markersize = 10, markeredgewidth=1.8, label="knots")

        ax.legend(fontsize='x-small',bbox_to_anchor=(0.91, 1), loc=2, borderaxespad=-1.)

    def __call__(self, *args):
        return self.pnt3D(args[0])

    def __repr__(self):
        return 'Nurbs curve:\n  degree: %s\n  cntrl: %s\n  uknots: %s' % (`self.degree`,`self.cntrl`,`self.uknots`)

class Line(Crv):
    """A straight line segment.

    >>> Line(p1=[0,0],p2=[1,1])
    Nurbs curve:
      degree: 1
      cntrl: array([[ 0.,  1.],
           [ 0.,  1.],
           [ 0.,  0.],
           [ 1.,  1.]])
      uknots: array([ 0.,  0.,  1.,  1.])
    """
    def __init__(self, p1 = (0,0,0), p2 = (1,0,0)):
        super(Line, self).__init__(np.transpose([p1,p2]), [0,0,1,1])
                       
class Polyline(Crv):
    """A polyline.

    >>> pnts = [[0,0],[5,2],[10,8]]
    >>> Polyline(pnts)
    Nurbs curve:
      degree: 1
      cntrl: array([[  0.,   5.,   5.,  10.],
           [  0.,   2.,   2.,   8.],
           [  0.,   0.,   0.,   0.],
           [  1.,   1.,   1.,   1.]])
      uknots: array([ 0. ,  0. ,  0.5,  0.5,  1. ,  1. ])
    """
    def __init__(self, pnts):
        pnts = np.transpose(np.asarray(pnts, np.float))
        npnts = pnts.shape[1]
        if npnts < 3:
            raise NURBSError, 'Point sequence error'
        cntrl = np.zeros((pnts.shape[0], 2 * npnts - 2), np.float)
        cntrl[:,0] = pnts[:,0]
        cntrl[:,-1] = pnts[:,-1]
        cntrl[:,1:-2:2] = pnts[:,1:-1]
        cntrl[:,2:-1:2] = pnts[:,1:-1]
        uknots = np.zeros(npnts * 2, np.float)
        uknots[0::2] = np.arange(npnts)
        uknots[1::2] = np.arange(npnts)
        super(Polyline, self).__init__(cntrl, uknots)

class UnitCircle(Crv):
    """NURBS representation of a unit circle in the xy plane.

    >>> c = UnitCircle()
    >>> assert c.degree==2
    >>> assert np.array_equiv(c.cntrl[3,1::2]/np.sqrt(2.)*2., np.ones(4))
    """
    def __init__(self):
        r22 = np.sqrt(2.)/2.
        uknots = [0., 0., 0., 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1., 1.,1.]
        cntrl = [[ 0.,  r22, 1., r22, 0., -r22, -1., -r22,  0.],
                 [-1., -r22, 0., r22, 1.,  r22,  0., -r22, -1.],
                 [ 0.,  0.,  0., 0.,  0.,  0.,   0.,  0.,   0.],
                 [ 1.,  r22, 1., r22, 1.,  r22,  1.,  r22,  1.]]
        super(UnitCircle, self).__init__(cntrl, uknots)

class Circle(UnitCircle):
    """NURBS representation of a circle in the xy plane
    with given radius (default = .5) and optional center.

    >>> c = Circle(radius=0.7, center=[1,0])
    >>> c.cntrl[1]
    array([-0.7       , -0.49497475,  0.        ,  0.49497475,  0.7       ,
            0.49497475,  0.        , -0.49497475, -0.7       ])
    """
    def __init__(self, radius = .5, center = None):
        super(Circle, self).__init__()
        if radius != 1.:
            self.trans(scale([radius, radius]))
        if center:
            self.trans(translate(center))

class Arc(Crv):
    """NURBS representation of an arc in the xy plane
    with given radius (default = 1.) and optional center,
    start angle (default = 0) and end angle. (default = 2*pi)

    >>> c0 = UnitCircle()
    >>> c1 = Arc()
    >>> u0 = np.linspace(0.25,0.5,num=5)
    >>> u1 = np.linspace(0.0,0.25,num=5)
    >>> np.allclose(c0.pnt3D(u0), c1.pnt3D(u1))
    True
    """
    def __init__(self, radius = 1.,center = None, sang = 0., eang = 2*math.pi):
        sweep = eang - sang # sweep angle of arc
        if sweep < 0.:
            sweep = 2.*math.pi + sweep
        if abs(sweep) <= math.pi/2.:
            narcs = 1   # number of arc segments
            knots = [0., 0., 0., 1., 1., 1.]
        elif abs(sweep) <= math.pi:
            narcs = 2
            knots = [0., 0., 0., 0.5, 0.5, 1., 1., 1.]
        elif abs(sweep) <= 3.*math.pi/2.:
            narcs = 3
            knots = [0., 0., 0., 1./3., 1./3., 2./3., 2./3., 1., 1., 1.]
        else:
            narcs = 4
            knots = [0., 0., 0., 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1., 1., 1.]

        dsweep = sweep/(2.*narcs);     # arc segment sweep angle/2
        
        # determine middle control point and weight
        wm = math.cos(dsweep)
        x  = radius*math.cos(dsweep)
        y  = radius*math.sin(dsweep)
        xm = x+y*math.tan(dsweep)

        # arc segment control points
        ctrlpt = np.array([[x, wm*xm, x], [-y, 0., y], [0., 0., 0.], [1., wm, 1.]], np.float)
        # build up complete arc from rotated segments
        coefs = np.zeros((4, 2*narcs+1), np.float)   # nurb control points of arc
        # rotate to start angle
        coefs[:,0:3] = np.dot(rotz(sang + dsweep), ctrlpt)
        xx = rotz(2*dsweep)
        for ms in range(2, 2*narcs,2):
            coefs[:,ms:ms+3] = np.dot(xx, coefs[:,ms-2:ms+1])
        if center:
            xx = translate(center)
            coefs = np.dot(xx, coefs)
        super(Arc, self).__init__(coefs, knots)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
