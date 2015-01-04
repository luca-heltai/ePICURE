import bpy
import sys
from mathutils import Vector
import os


def clear_scene():

	if(len(bpy.context.selected_objects) != 0):
		bpy.ops.object.select_all(action="TOGGLE")

	bpy.ops.object.select_all(action="TOGGLE")

	bpy.ops.object.delete()

def MakePolyLine(objname, curvename, cList):
	curvedata = bpy.data.curves.new(name=curvename, type='CURVE')
	curvedata.dimensions = '3D'

	objectdata = bpy.data.objects.new(objname, curvedata)
	objectdata.location = (0,0,0) #object origin
	bpy.context.scene.objects.link(objectdata)

	polyline = curvedata.splines.new('NURBS')
	polyline.points.add(len(cList)-1)
	for num in range(len(cList)):
		polyline.points[num].co = (cList[num])+(w,)

	polyline.order_u = len(polyline.points)-1
	polyline.use_endpoint_u = True


bpy.ops.wm.addon_enable(module="game_engine_save_as_runtime")
bpy.ops.wm.save_as_runtime(player_path=os.environ['BLENDER_PLAYER'], filepath=os.environ['SAVEPATH']+os.environ['FILENAME'], copy_python = True, overwrite_lib = False)




"""
class CompositeObject(object):

	def __init__(self, m, n, ntimeframes):
	#CompositeObject, initialized with m, n, ntimeframes or directly with two filenames

	@classmethod
    def fromfilename(cls, file_shapes, file_positions):



		return cls(m,n,ntimeframes)

	def read_shapes(filename):

	def read_positions(filename):

	def pose_object(shape, position):

		listOfVectors = [(0,0,0),(1,0,0),(2,0,0),(2,3,0),(0,2,1)]

		MakePolyLine("NameOfMyCurveObject", "NameOfMyCurve", listOfVectors)

"""
