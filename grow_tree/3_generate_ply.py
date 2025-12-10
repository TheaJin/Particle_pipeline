import os
import sys
import pymesh
import annotate as dp

data_path = '/hpc/gjin250/mapclient/output'
study = 'virtual_child'
subject =  '025'
#volume = 'Pre'

subject_path = os.path.join(data_path, study, subject)
geometry_path = os.path.join(subject_path)#, 'Lobe', 'Lobe_mesh')

LUL_stl = os.path.join(geometry_path, 'LUL' + '.stl')
LLL_stl = os.path.join(geometry_path, 'LLL' + '.stl')
RUL_stl = os.path.join(geometry_path, 'RUL' + '.stl')
RML_stl = os.path.join(geometry_path, 'RML' + '.stl')
RLL_stl = os.path.join(geometry_path, 'RLL' + '.stl')

LUL_mesh = pymesh.load_mesh(LUL_stl)
LLL_mesh = pymesh.load_mesh(LLL_stl)
RUL_mesh = pymesh.load_mesh(RUL_stl)
RML_mesh = pymesh.load_mesh(RML_stl)
RLL_mesh = pymesh.load_mesh(RLL_stl)


resolution = "low"

LUL_surf, LUL_vert, target_LUL = dp.py_mesh_and_repair(LUL_mesh, resolution) # returns unscaled data
LLL_surf, LLL_vert, target_LLL = dp.py_mesh_and_repair(LLL_mesh, resolution) # returns unscaled data
RUL_surf, RUL_vert, target_RUL = dp.py_mesh_and_repair(RUL_mesh, resolution) # returns unscaled data
RML_surf, RML_vert, target_RML = dp.py_mesh_and_repair(RML_mesh, resolution) # returns unscaled data
RLL_surf, RLL_vert, target_RLL = dp.py_mesh_and_repair(RLL_mesh, resolution) # returns unscaled data

# Translate nodes (Optional)
translation = False #True
if translation:
    translation_vector = [-169.7, -306.2, 937.7] # For alfred3 from Metadata
    LUL_surf.points += translation_vector
    LLL_surf.points += translation_vector
    RUL_surf.points += translation_vector
    RML_surf.points += translation_vector
    RLL_surf.points += translation_vector



### Export poly mesh files ###

LLL_surf.save(os.path.join(geometry_path, 'LLL_surf.ply'), binary=False)
LUL_surf.save(os.path.join(geometry_path, 'LUL_surf.ply'), binary=False)
RML_surf.save(os.path.join(geometry_path, 'RML_surf.ply'), binary=False)
RLL_surf.save(os.path.join(geometry_path, 'RLL_surf.ply'), binary=False)
RUL_surf.save(os.path.join(geometry_path, 'RUL_surf.ply'), binary=False)