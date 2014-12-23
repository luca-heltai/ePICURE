import numpy as np
from utilities.least_square import *

def test_least_square():
    vector_space = UniformLagrangeVectorSpace(4)
    function = lambda x: np.sin(10*np.pi*x)
    points = np.array([0, 0.3, 1.0])
    works = True
    try:
        least_square(function, vector_space, points)
        works = False
    except:
        pass
    assert works, \
           'Expecting failure!'

    points = np.array([0, 0.3, 0.6, 1.0])
    works = False
    try:
        least_square(function, vector_space, points)
        works = True
    except:
        pass
    assert works, \
           'Unexpected failure! The number of points is equal to the dofs'

    points = np.array([0, 0.3, 0.6, 0.8, 1.0])
    works = False
    try:
        least_square(function, vector_space, points)
        works = True
    except:
        pass
    assert works, \
           'Unexpected failure! The number of points is greater than the dofs'

def vectorial_test_least_square():
    vector_space = UniformLagrangeVectorSpace(4)
    function=lambda t: np.array([t,t**2,t**3]).T
    points = np.array([0, 0.3, 0.6, 0.8, 1.0])
    works = False
    try:
        least_square(function, vector_space, points)
        works = True
    except:
        pass
    assert works, \
           'Unexpected failure in the vector case!'

    n=2**3+1
    knots = np.linspace(0,1,n)
    vector_space = IteratedVectorSpace(UniformLagrangeVectorSpace(5),knots)
    function=lambda t: np.array([t,t**2,t**3]).T
    points = np.linspace(0,1,2**7+1)
    eps=10e-1
    assert (abs(least_square(function, vector_space, points)[17][:]-np.array([.5,.25,.125]))<eps).all(),\
            'Unexpected failure! Bad approximation'
