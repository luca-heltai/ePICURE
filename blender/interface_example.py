import sys
import os
import bpy
sys.path.insert(0, os.environ['CURRENTDIR'])

from blender_interface import *

#Open sample blender file
open_blender_file("tracking.blend")
#clear_scene()
#Enable AnimALL module
enable_AnimALL()

#Load tracked coordinates from the file and (just for this simple example) save the 3D coordinates 
#of the points in coords.
# coords = sample_ortogonal_camera_solver(get_tracked_coords())

#Save the coords for each frame inside a file
# write_to_file_coords(coords,filename="sample.shape")

#Create a composite object
# swimmer = CompositeObject.load_from_file("sample.shape","sample.pos","swimmer")

#Place the object as in the frame 0 in the scene.
#swimmer.pose_object()

#Animate the object.
#animate_CO(swimmer)
