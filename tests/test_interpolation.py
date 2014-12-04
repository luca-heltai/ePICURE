import numpy as np
from utilities.interpolation import *

def test_interpolation():
    vector_space = UniformLagrangeVectorSpace(3)
    function = lambda x: np.sin(10*np.pi*x)
    points = np.array([0, 0.3, 1.0])
    c = interpolation(function, vector_space, points)
    works = True
    for i in range(3):
        s = 0
        for j in range(3):
            s += c[j]*vector_space.basis(j)(points[i])
        if abs(s - function(points[i])) > 1e-15:
            works = False
            break
    assert works, \
           'The interpolation and the function assume different values on the interpolation points!'

    points = np.array([0, 0.3, 0.6, 1.0])
    works = True
    try:
        interpolation(function, vector_space, points)
        works = False
    except:
        pass
    assert works, \
           'Expecting failure when the number of points is greater than the dofs!'

    points = np.array([0, 1.0])
    works = True
    try:
        interpolation(function, vector_space, points)
        works = False
    except:
        pass
    assert works, \
           'Expecting failure when the number of points is less than the dofs!'
