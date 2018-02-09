import numpy as np
from numpy import linalg as LA
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from utilities.matrices import *
import interfaces
from utilities import least_square
from .arclength import *
from scipy.interpolate import interp1d

class ArcLenghtCurve(object):
    """A python utility used to construct and manipulate an arclength curve
    in the space."""

    def __init__(self,vs):
        self.vs = vs
        self.space_dim = 0 #to be implemented
        self.gamma = ''
        self.coords = ''

    def from_lambda_to_coords(self, s_space=np.linspace(0,1,1025)):
        """This method allows us to pass from lambda representation to
        coords representation of self.gamma."""
        assert (self.gamma.size != 0), \
                "self.gamma is not defined."
        # calculate the dimension of the space given the curve
        gamma = lambda s: self.gamma(s).T
        np.shape(gamma(self.s_space))
        self.coords = least_square(gamma, self.vs, s_space)

    def from_coords_to_lambda(self):
        """This method allows us to pass from coords representation to
        lambda representation of self.gamma."""
        assert (self.coords.size != 0), \
               "self.coords is not defined."
        self.gamma = self.vs.element(self.coords)

    def first_derivatives(self):
        """ If self.gamma is defined this method returns first, second,
        and third derivative of the curve."""
        assert (self.coords.size != 0), \
            "self.coords is not defined."
        self.ds   = self.vs.element_der(self.coords,1)
        self.dds  = self.vs.element_der(self.coords,2)
        self.ddds = self.vs.element_der(self.coords,3)

    def first_derivative_modulus(self):
        assert (self.gamma.size != 0), \
            "self.gamma is not defined."
        self.first_derivatives()
        print (((self.ds(self.s_space).T).dot(self.ds(self.s_space))).diagonal())

    def first_normal_modulus(self):
        print (((self.normal(self.s_space)).dot(self.normal(self.s_space).T)).diagonal())

    def torsion(self):
        """This method provides a lambda function representing the torsion
            of self.gamma."""
        self.first_derivatives()
        def tau(t):
            DX = self.ds(t)
            DDX = self.dds(t)
            DDDX = self.ddds(t)
            den = np.array(LA.norm(np.cross(DX.T,DDX.T).T,axis=0)**2)
            ids = np.array(den==0)
            den[ids]+=np.finfo(np.float64).eps
            return np.sum(np.cross(DX.T,DDX.T).T*DDDX,axis=0)/den
        self.tau = lambda t: tau(t)

    def curvature(self):
        """This method provides a lambda function representing the curvature
            of self.gamma."""
        self.first_derivatives()
        def curv(t):
            DX = self.ds(t)
            DDX = self.dds(t)
            den = np.array((LA.norm(DX, axis=0))**3)
            ids = np.array(den==0)
            den[ids]+=np.finfo(np.float64).eps
            return LA.norm(np.cross(DX.T,DDX.T).T,axis=0)/den
        self.kappa = lambda t: curv(t)

    def reparametrise(self):
        """This method reparametrise the curve in order to find an
        arclenght curve."""
        if(self.coords==''):
            self.from_lambda_to_coords()
        reparamCurve = ArcLengthParametrizer(self.vs,self.coords,\
                            arcfactor = int(len(self.s_space))*100)
        self.gamma = lambda t: reparamCurve.curve(t)

    def curve_from_curvature(self,\
                    start = 0, end = 1, \
                    x0 = np.array([0,0,0]),\
                    frenet0 = np.eye(3,3).reshape((-1,)),\
                    L = 1):
        """This method provides a curve starting from curvature and torsion."""
        k = lambda s : np.array([self.kappa(s)])
        t = lambda s : np.array([self.tau(s)])

        assert np.shape(k(0)) == np.shape(t(0)), \
                "ERROR: Curvature and Torsion has different shape as numpy array."

        try:
        # np.max is needed to avoid error using column/row array.
            dim = np.max(np.shape(k(0)))
        except:
            dim = 1

        KAPPA = np.array([[0,1,0],[-1,0,0],[0,0,0]])
        TAU = np.array([[0,0,0],[0,0,1],[0,-1,0]])
        K = lambda s : (KAPPA.reshape(-1,1)*k(s)).reshape(-1,)
        T = lambda s : (TAU.reshape(-1,1)*t(s)).reshape(-1,)
        KT = lambda s : (K(s)+T(s)).reshape(3,3,-1).T.transpose(0,2,1)

        # Default value for frenet initial setting:
        # It is the identity matrix for every curve
        onesk = np.array(np.ones(dim))
        frenet0 = (np.eye(3).reshape(-1,1)*onesk).reshape(3,3,-1).T.transpose(0,2,1).reshape(dim,-1).T.reshape(-1,)

        def Matrix(frenet,s):
            MMR = lambda s: [KT(s)[i].dot(frenet.reshape(3,3,dim).T.transpose(0,2,1)[i]) for i in range(dim)]
            #return np.array(zip(MMR(s))).reshape(dim,-1).T.reshape(-1,)
            return np.array(MMR(s)).reshape(dim,-1).T.reshape(-1,)

        frenet = odeint(Matrix, frenet0, self.s_space)

        F = frenet.reshape(-1,3,3,dim).transpose(0,2,1,3)

        # Decomposing Frenet-Serret set
        L , N , B = F[:,:,:,:].transpose(2,3,1,0)
        curve_discrete = np.cumsum(np.squeeze(L),axis=-1)/len(self.s_space)
        self.gamma = lambda t: np.array(interp1d(self.s_space,curve_discrete)(t))
        self.normal = lambda t: np.array(interp1d(self.s_space,np.squeeze(N))(t))
        self.binormal = lambda t: np.array(interp1d(self.s_space,np.squeeze(B))(t))
        self.reparametrise()

    def plot(self, \
            binormal=False,  	binormal_lenght=.0005, \
            normal=False, 		normal_lenght=.0005):
        """This method plot self.gamm in a 3D plot. Options that an be used
        to analyse this plot are the following:
        - normal = False/True : this option allows to plot the normal
            vector fiels;
        - normal_lenght : it is the lenght of each arrow of the vector field;
        - binormal = False/True : this option allows to plot the binormal
            vector fiels;
        - binormal_lenght : it is the lenght of each arrow of the vector
            field."""
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        plt.axis('equal')
        ax.plot(self.gamma(self.s_space)[0], self.gamma(self.s_space)[1], self.gamma(self.s_space)[2],'r', label='curve')
        if normal:
            ax.quiver(\
                self.gamma(self.s_space)[0]+normal_lenght*self.normal(self.s_space)[0], \
                self.gamma(self.s_space)[1]+normal_lenght*self.normal(self.s_space)[1], \
                self.gamma(self.s_space)[2]+normal_lenght*self.normal(self.s_space)[2], \
                +self.normal(self.s_space)[0], +self.normal(self.s_space)[1], +self.normal(self.s_space)[2], \
                length=normal_lenght,  color='blue')
        if binormal:
            ax.quiver(\
                self.gamma(self.s_space)[0]+binormal_lenght*self.binormal(self.s_space)[0], \
                self.gamma(self.s_space)[1]+binormal_lenght*self.binormal(self.s_space)[1], \
                self.gamma(self.s_space)[2]+binormal_lenght*self.binormal(self.s_space)[2], \
                +self.binormal(self.s_space)[0], +self.binormal(self.s_space)[1], +self.binormal(self.s_space)[2], \
                length=binormal_lenght,  color='green')
        ax.legend()
        plt.show()
        # plt.close()

class ALCFromCoords(ArcLenghtCurve):
    """	Given a set of coordinates with respect to the basis
        of the vector space, this methot returns a lambda function
        representing the curve in the space."""
    def __init__(self, vs, coords):
        self.vs = vs
        self.coords = coords
        self.s_space = np.linspace(0.0,1.0,num=1025)
        self.from_coords_to_lambda()

class ALCFromLambda(ArcLenghtCurve):
    """ Given a vectorial lambda function this method returns the coordinates
        with respect to the basis of the vector space."""
    def __init__(self, vs, gamma, s_space=np.linspace(0.0,1.0,num=1025)):
        self.coords = ''
        self.vs = vs
        self.s_space=s_space
        self.gamma = lambda t : gamma(t)
        self.from_lambda_to_coords()
        self.curvature()
        self.torsion()

class ALCFromKappaAndTau(ArcLenghtCurve):
    """ Given two lambda functions representing curvature and torsion,
        this class returns an arclenght curve. Moreover, curvature
        and torsione are recalculated with respect the new curve."""
    def __init__(self, vs, kappa, tau,\
                    x0 = np.array([0,0,0]),\
                    frenet0 = np.eye(3,3).reshape((-1,)),\
                    L = 1,\
                    s_space = np.linspace(0.0,1.0,num=1025)):
        self.coords = ''
        self.kappa = kappa
        self.tau = tau
        self.vs = vs
        self.x0 = x0
        self.frenet0 = frenet0
        self.L = L
        self.s_space=s_space
        self.curve_from_curvature(
                    start = 0, end = 1, \
                    x0 = self.x0,\
                    frenet0 = self.frenet0,\
                    L = self.L)
        self.from_lambda_to_coords()
        self.curvature()
        self.torsion()