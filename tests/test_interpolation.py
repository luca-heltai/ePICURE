from utilities.interpolation import *

def test_interpolation():
    vector_space = UniformLagrangeVectorSpace(3)
    function = lambda x: sin(10*pi*x)
    points = np.array([0, 0.3, 1.0])
    c = interpolation(function, vector_space, points)
    works = True
    for i in range(3):
        s = 0
        for j in range(3):
            s += c[j]*vector_space.basis(j)(points[i])
        if abs(s - funz(points[i])) > 1e-15:
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
           'Expecting failure!'

    points = np.array([0, 1.0])
    works = True
    try:
        interpolation(function, vector_space, points)
        works = False
    except:
        pass
    assert works, \
           'Expecting failure!'
