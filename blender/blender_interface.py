import bpy
import sys
import os
import numpy as np

"""
This is a blender interface able to import inside blender a NURBS object defined from two input files.
This interface is also able to save the coordinates for each frame of tracked tracks using blender tools.
"""
__author__ =  'Leonardo Romor'
__version__=  '1.0'



################################################################################
#	BLENDER UTILS
################################################################################



def open_blender_file(file):
	"""
	Simple function to load a blender file.

	input:
		file - (String)the filepath with the filename.
	"""

	bpy.ops.wm.open_mainfile(filepath=file)



def Saveruntime():
	"""
	Function that automatically create a standalone application
	Requires the following path variables to be defined (for example using export BLENDER_PLAYER=/..):

		- BLENDER_PLAYER = the path to the blender standalone player
		- SAVE_PATH = the path in which you wish to save the file
		- FILENAME = the name of the file
	"""

	#PATHS
	blender_player=os.environ['BLENDER_PLAYER']
	save_path=os.environ['SAVEPATH']
	filename=os.environ['FILENAME']

	bpy.ops.wm.addon_enable(module="game_engine_save_as_runtime")
	bpy.ops.wm.save_as_runtime(player_path=blender_player, filepath=save_path+filename, copy_python = True, overwrite_lib = False)



def enable_AnimALL():
	"""
	Enable the blender module AnimALL that allow us to keyframe even vertices
	"""

	bpy.ops.wm.addon_enable(module="animation_animall")



################################################################################
#	MOTION TRACKING
################################################################################



def get_tracked_coords():
	"""
	From the loaded blender file, if they exist, this function return a list containing all the coordinates
	of the tracks for each frame and for each movieclip.

	example:
		coords[(#clip)][#track][0] = takes the x coordinate from n#clip of the n#track
		coords[(#clip)][#track][1] = takes the y coordinate from n#clip of the n#track

	usually the tracked coordinates range from 0 to 1. (0,0) represent the bottom-left angle.
	"""

	coords=[]
	#Loop in movievlips
	for idxclip,clip in enumerate(bpy.data.movieclips):
		coords.append([])
		#Loop in tracks
		for idxtrack,track in enumerate(clip.tracking.tracks):
			coords[idxclip].append([])
			frameno = 0
			while True:
				markerAtFrame = track.markers.find_frame(frameno)
				if not markerAtFrame:
					break
				frameno += 1
				coords[idxclip][idxtrack].append([markerAtFrame.co[0],markerAtFrame.co[1]])

	return coords



def write_to_file_coords(coords,filename="coords.shape"):
	"""
	Save for the input coords from the simple camera solver into a .shape file which contains
	blocks separated by a newline. Each of the blocks represent one frame, every line of the block is one track.
	
	input:
		- coords = Variable which contains sample_ortogonal_camera_solver() output.
		- filename = output filename.
	"""

	out_file = open(filename,"w")
	for index,frame in enumerate(coords):
		for track in frame:
			out_file.write(str(track[0])+" "+str(track[1])+" "+str(track[2])+"\n")

		out_file.write("\n")



def sample_ortogonal_camera_solver(coords):
	"""
	(SAMPLE)
	A very simple camera solver which takes from two "ortoghonal" ortographics cameras the coordinates of n tracked objects
	for example camera1 display xz plane, camera2 yz.

	input:
		- coords = accept as input the output of get_tracked_coords() for this special case.
	"""

	frames=[]
	#X,Z
	camera1=coords[0]
	for i in range(0, len(camera1[0])-1):
		frame=[]
		for tracks in camera1:
			frame.append( [tracks[i][0]-0.5,0,tracks[i][1]-0.5] )

		frames.append(frame)		
	#Y
	camera2=coords[1]
	for i in range(0, len(camera2[0])-1):
		for idx,tracks in  enumerate(camera2):
			frames[i][idx][1]=tracks[i][0]-0.5

	##FRAME[TRACK[X,Y,Z]]
	return frames



################################################################################
#	SCENE MANAGEMENT
################################################################################



def clear_scene():
	"""
	A simple function to clear the scene.
	"""

	if(len(bpy.context.selected_objects) != 0):
		bpy.ops.object.select_all(action="TOGGLE")

	bpy.ops.object.select_all(action="TOGGLE")

	bpy.ops.object.delete()



#-------	OBJECT MANIPULATION



def animate_CO(CompositeObject):
	"""
	Animate the shape of the object. The position and rotation animation of the object is still not implemented.

	input:
		- CompositeObject = An instance of CompositeObject.
	"""

	blend_object = bpy.data.objects[CompositeObject.name]
	
	cu=blend_object.data

	points=cu.splines[0].points
	bpy.data.window_managers["WinMan"].key_points=True

	bpy.context.scene.objects.active = blend_object

	for index, frame in enumerate(CompositeObject.object_shape):
		if(index==0):
			bpy.ops.anim.insert_keyframe_animall()
		else:
			bpy.context.scene.frame_set(index)
			bpy.ops.anim.insert_keyframe_animall()
			for pidx,point in enumerate(points):
				bpy.ops.anim.insert_keyframe_animall()
				point.co[0]=frame[pidx][0]
				bpy.ops.anim.insert_keyframe_animall()
				point.co[1]=frame[pidx][1]
				bpy.ops.anim.insert_keyframe_animall()
				point.co[2]=frame[pidx][2]
	
	bpy.ops.wm.save_mainfile()



def MakePolyLine(objname, curvename, cList,cposList):
	"""
	A function that create and place a NURBS curve from a set of points and weights.

	input:
		- objname = name of the final object
		- curvename = name of the curve inside the object
		- cList = list containing all the points in sequence
		- w = the weight of the NURBS
	"""

	w=90 #this is chosen as default since I don't have yet a routing which evaluates the control points of the nurbs.
	curvedata = bpy.data.curves.new(name=curvename, type='CURVE')
	curvedata.dimensions = '3D'

	objectdata = bpy.data.objects.new(objname, curvedata)
	objectdata.location = (cposList[0][0],cposList[0][1],cposList[0][2]) #object origin
	bpy.context.scene.objects.link(objectdata)

	polyline = curvedata.splines.new('NURBS')
	polyline.points.add(len(cList)-1)
	for num in range(len(cList)):
		polyline.points[num].co = (cList[num])+[w]

	polyline.order_u = len(polyline.points)-1
	polyline.use_endpoint_u = True



class CompositeObject(object):
	"""
	In this version of the interface a CompositeObject can only be an object containing a NURBS curve.
	A CompositeObject can be inizialized in two different ways, from a file or directly with two lists.

	member variables:
		- object_position = describes the object position for each frame.
		- object_shape = describes the object shape for each frame.
		- name = the name of the object
	"""

	def __init__(self,object_shape,object_position,name):
		"""
		To instance a CompositeObject directly from two lists. 

		input:
			- object_shape = list shaped in this way: object_shape[frame][point][xyz]
			- object_position = list shaped in this way: object_position[frame][0][xyz]
		"""

		self.object_position = object_position
		self.object_shape = object_shape
		self.name = name


	@classmethod
	def load_from_file(cls,file_shapes,file_position,name):
		"""
		To instance a CompositeObject from a file.

		exampe:
			swimmer = CompositeObject.load_from_file("sample.shape","sample.pos","swimmer")
		input:
			- file_shapes = file containing blocks representing each frame, where each line represent the coordinate
				of a control vertex in that frame.
			- file_position = as for file_shape but this time each blocks has just one which represent the coordinate of the global
				position of the object.
		"""

		object_shape=CompositeObject.shape_file_parser(file_shapes)
		object_position=CompositeObject.position_file_parser(file_position)

		return cls(object_shape,object_position,name)

	
	@staticmethod	
	def shape_file_parser(filename):
		"""
		This function parse the shape file.

		input:
			- filename = path to the input file.
		"""

		#SPLITTING THE FILE IN BLOCKS
		empty_lines = 0
		blocks=[]
		for line in open(filename):
			# Check for empty/commented lines
			if line.startswith('#') or line.startswith('\n'):
	        	# If 1st one: new block
				if empty_lines == 0:
					blocks.append([])
				empty_lines += 1

	    	# Non empty line: add line in current(last) block
			else:
				empty_lines = 0
				if(len(blocks)==0):
					blocks.append([])
				blocks[-1].append(line)
		if(len(blocks[-1])==0):
			blocks.pop()

		#SAVING ALL THE BLOCKS(FRAME) SUCH THAT object_shape[FRAME][POINT][XYZ]
		ntimeframes=len(blocks)
		object_shape=[]
		for indx,block in enumerate(blocks):
			object_shape.append([])
			for line in block:
				array = [float(x) for x in line.split()] 
				object_shape[indx].append(array)

		return object_shape


	@staticmethod	
	def position_file_parser(filename):
		"""
		This function parse the position file.

		input:
			- filename = path to the input file.
		"""

		empty_lines = 0
		blocks=[]
		for line in open(filename):
			# Check for empty/commented lines
			if line.startswith('#') or line.startswith('\n'):
	        	# If 1st one: new block
				if empty_lines == 0:
					blocks.append([])
				empty_lines += 1

	    	# Non empty line: add line in current(last) block
			else:
				empty_lines = 0
				if(len(blocks)==0):
					blocks.append([])
				blocks[-1].append(line)
		if(len(blocks[-1])==0):
			blocks.pop()

		ntimeframes=len(blocks)
		object_shape=[]
		for indx,block in enumerate(blocks):
			object_shape.append([])
			for line in block:
				array = [float(x) for x in line.split()] 
				object_shape[indx].append(array)

		return object_shape


	def pose_object(self,frame=0):
		"""
		This function place the object in the blender scene in the current frame from the specified "(shape,position)" frame.

		input:
			- frame = the Object "internal" frame to place in the current blender frame.
		"""
		
		MakePolyLine(self.name, self.name ,self.object_shape[frame],self.object_position[frame])

