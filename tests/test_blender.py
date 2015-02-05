import sys

sys.path.append('../utilities')
sys.path.append('../interfaces')
sys.path.append('../blender')
from blender_interface import *
from nose.tools import *
import numpy as np 



def test_composite_object():

	shape=[[],[]]
	shape[0]=[0,232,241.0]
	shape[1]=[0,12,152.0]

	pos=[[],[]]
	pos[0]=[0,232,241.0]
	works = False
	try:
		swimmer = CompositeObject(shape,pos,"swimmer")
		works=True
	except:
		pass
	assert works

