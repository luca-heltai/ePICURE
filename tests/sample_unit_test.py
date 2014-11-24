import numpy as np

def test_point_in_sequence():
    # make sure a point is in the sequence
    p = .12
    x = np.linspace(0,1,11)
    y = x[(x<p)]
    assert len(y) ==  2
