#!/usr/bin/env python

import argparse
import os
import shutil
#import matplotlib.pyplot as plt
import numpy as np

from aether.diagnostics import set_diagnostics_on
from aether.exports import export_1d_elem_geometry, export_node_geometry, export_elem_field, export_1d_elem_field, \
    export_terminal_solution,export_node_field
from aether.imports import import_ventilation
from aether.indices import define_problem_type, get_ne_radius
from aether.geometry import add_mesh, append_units, define_node_geometry, define_1d_elements, define_rad_from_file, \
    define_rad_from_geom, refine_1d_elements, reorder_tree, set_initial_volume
from aether.geometry import list_tree, scale_tree, scale_radii
from aether.particle_transport import solve_particles_decoupled
from aether.ventilation import evaluate_uniform_flow
from aether.field_utilities import scale_flow_to_inlet

# def write_parameters(subject):
#     # dictionary to contain the subject lung measurement data
#     subject_data = { 'P2BRP001-H653': ['f', 20, 1.63, 54, 20.45, 5.64, 2.8],
#                      'P2BRP021-H673': ['f', 23, 1.8, 70.45, 21.74, 7.44, 3.38],
#                      'P2BRP030-H682': ['f', 21, 1.7, 62.1, 21.49, 5.75, 3.6],
#                      'P2BRP032-H684': ['f',27, 1.63, 64.5, 24.28, 5.08, 2.77],
#                      'P2BRP042-H18': ['f', 24, 1.6, 56, 21.875, 4.06, 1.77],
#                      'P2BRP076-H1335': ['f', 20, 1.73, 66.8, 22.35, 6.57, 3.45],
#                      'P2BRP115-H5977': ['m', 21, 1.73, 85.2, 28.43, 7.41, 2.45],
#                      'P2BRP139-H6229': ['f', 23, 1.568, 52.6, 21.39, 3.98, 1.66],
#                      'P2BRP143-H6310': ['f', 20, 1.56, 52, 21.37, 4.48, 2.06],
#                      'P2BRP242-H11303': ['m', 22, 1.821, 82.6, 24.91, 8.35, 4.01],
#                      'P2BRP243-H11750': ['m', 21, 1.78, 80.7, 25.53, 7.09, 3.62],
#                      'P2BRP247-H11779': ['m', 24, 1.78, 83.1, 26.23, 8.86, 4.38],
#                      'P2BRP268-H12816': ['m', 23, 1.87, 80.9, 23.11, 7.13, 3.4] }
#     
#     opf = open('Parameters/params_evaluate_flow.txt', "w+")
#     opf.write("FRC        %4.2f\n" % (subject_data[subject][6]) )
#     opf.write("constrict  1.0\n")
#     opf.write("T_interval 4.0\n")
#     opf.write("Gdirn      3\n") 
#     opf.write("press_in   0.0 \n")
#     opf.write("COV        0.1\n") 
#     opf.write("RMaxMean   1.29\n") 
#     opf.write("RMinMean   0.79\n") 
#     opf.write("i_to_e_ratio 1.0\n") 
#     opf.write("refvol     0.5\n") 
#     opf.write("volume_target 500000\n") 
#     opf.write("pmus_step  -196.133\n") 
#     opf.write("expiration_type active \n")
#     opf.write("chest_wall_compliance 2039.4324\n") 
#     opf.close
# 
#     return subject_data[subject][6] * 1.0e+6

def HLA_subject_data():
    # dictionary to contain the subject lung measurement data
    # for H1335, Rtrach = hght * 0.74394657607831 + 1.63* 3.1634738580296 = 6.44348996520372
    #                               sex   age   hght   wght    BMI    TLC     FRC   Rtrach
    subject_data = {'P2BRP001-H653': ['f', 20, 1.63, 54, 20.45, 5.64, 2.8, 6.4571],
                    'P2BRP021-H673': ['f', 23, 1.8, 70.45, 21.74, 7.44, 3.38, 6.4571],
                    'P2BRP030-H682': ['f', 21, 1.7, 62.1, 21.49, 5.75, 3.6, 5.7669],
                    'P2BRP032-H684': ['f', 27, 1.63, 64.5, 24.28, 5.08, 2.77, 6.5260],
                    'P2BRP042-H18': ['f', 24, 1.6, 56, 21.875, 4.06, 1.77, 0.0],
                    'P2BRP076-H1335': ['f', 20, 1.73, 66.8, 22.35, 6.57, 3.45, 6.4435],
                    'P2BRP115-H5977': ['m', 21, 1.73, 85.2, 28.43, 7.41, 2.45, 7.9411],
                    'P2BRP139-H6229': ['f', 23, 1.568, 52.6, 21.39, 3.98, 1.66, 5.8750],
                    'P2BRP143-H6310': ['f', 20, 1.56, 52, 21.37, 4.48, 2.06, 5.2187],
                    'P2BRP242-H11303': ['m', 22, 1.821, 82.6, 24.91, 8.35, 4.01, 7.7962],
                    'P2BRP243-H11750': ['m', 21, 1.78, 80.7, 25.53, 7.09, 3.62, 7.3996],
                    'P2BRP247-H11779': ['m', 24, 1.78, 83.1, 26.23, 8.86, 4.38, 6.9213],
                    'P2BRP268-H12816': ['m', 23, 1.87, 80.9, 23.11, 7.13, 3.4, 7.7227]}
    return subject_data

# -----------------------------------------------------------------------------------------------------

def write_parameters(subject_data, subject):
    opf = open('Parameters/params_evaluate_flow.txt', "w+")
    opf.write("FRC        %4.2f\n" % (subject_data[subject][6]))
    opf.write("constrict  1.0\n")
    opf.write("T_interval 4.0\n")
    opf.write("Gdirn      3\n")
    opf.write("press_in   0.0 \n")
    opf.write("COV        0.1\n")

    # zeroG
    # opf.write("RMaxMean   1.0\n")
    # opf.write("RMinMean   1.0\n")

    # oneG upright
    opf.write("RMaxMean   1.29\n")
    opf.write("RMinMean   0.79\n")

    # oneG invert
    # opf.write("RMaxMean   0.79\n")
    # opf.write("RMinMean   1.29\n")

    opf.write("i_to_e_ratio 1.0\n")
    opf.write("refvol     0.5\n")
    opf.write("volume_target 500000\n")
    opf.write("pmus_step  -196.133\n")
    opf.write("expiration_type active \n")
    opf.write("chest_wall_compliance 2039.4324\n")
    opf.close

    # scaling from FRC to TLC includes 1 L of tissue (FRC and TLC are air)
    scale_to_FRC = (subject_data[subject][6] + 1.0) / (subject_data[subject][5] + 1.0)
    return subject_data[subject][7], scale_to_FRC


def main():
    parser = argparse.ArgumentParser(description='Run particle transport for a selected subject')
    parser.add_argument('-s','--subject', help='Input the subject name',required=True)
    parser.add_argument('-p','--particle_size', help='Input the particle size in microns',required=True)
    parser.add_argument('-g','--gravity', help='Input the gravity direction (0/1/2/3/-3)',required=True)
    parser.add_argument('-v','--vent', help='Input the ventilation file code',required=True)
    
    args = parser.parse_args()
    # new
    subject_data = HLA_subject_data()
    abbrv_subject = args.subject
    for key in subject_data:
        if abbrv_subject in key:
            subject = key
    
    #subject = args.subject
    print(subject)
    particle_size = float(args.particle_size)
    Gdirn = int(args.gravity)
    vent = args.vent
    volume = 'FRC'

    set_diagnostics_on(False)

    export_directory = os.path.join('output', str(abbrv_subject))
    if not os.path.exists(export_directory):
        os.makedirs(export_directory)
    plot_me = False

    define_problem_type('particle_transport')
    trachea_rad, scale_to_FRC = write_parameters(subject_data, subject)
    filename = subject + '_grown_airways'

    define_node_geometry('geometries/' + filename)
    define_1d_elements('geometries/' + filename)


    # scale the TLC airway geometries to approximate FRC lengths
    scale_f = scale_to_FRC ** (1.0 / 3.0)
    scale_tree('x', scale_f)
    scale_tree('y', scale_f)
    scale_tree('z', scale_f)

    # --reading from file based on CT: giving airways that are too narrow!
    
    # TJ - param to scale the deadspace size, 1.0 = unscaled
    define_rad_from_file('geometries/upper_template.ipfiel', 1.0)
    scale_trachea = trachea_rad / 8.0  # divide by the trachea size in the template file
    scale_radii(scale_trachea)  # scale all of the template by the same value

    min_rad = 0.22
    define_rad_from_geom('fit', 0.0, 'inlet', min_rad)
    list_tree()

    refine_elements = list(range(1,21))
    refine_1d_elements(refine_elements, 10)
    reorder_tree()

    append_units()

    COV = 0.1
    if Gdirn > 0:
        Rmax = 1.29
        Rmin = 0.79
    elif Gdirn < 0:
        Rmax = 0.79
        Rmin = 1.29
    else:
        Rmax = 1.0
        Rmin = 1.0


    # Set the working directory to the this files directory and then reset after running simulation.
    file_location = os.path.dirname(os.path.abspath(__file__))
    cur_dir = os.getcwd()
    os.chdir(file_location)

    # read in results of ventilation model
    import_filename = export_directory + '/' + subject + '_ventilation_field' # + vent
    import_ventilation(import_filename + '.exelem')
    
    
    # this is a normalised field so needs to be scaled to our required flow value
    tidal_volume = 500000  # mm3
    time_inspiration = 2.0 # seconds
    flow_rate = tidal_volume / time_inspiration

    # scale flow to alveolar ventilation
    scale_flow_to_inlet(flow_rate,'V')

    initial_concentration = 0
    inlet_concentration = 0.001  # consistent with inner workings of lungsim

    # solve particle transport
    solve_particles_decoupled(initial_concentration, inlet_concentration, particle_size)

    name = "{0}_p{1}".format(subject, int(particle_size))
    #filename = export_directory + '/' + name
    filename = os.path.join(export_directory, str(particle_size), name)
    groupname = 'vent_model'

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    export_node_geometry(filename + '.exnode', groupname)
    export_1d_elem_geometry(filename + '.exelem', groupname)

    field_name = 'lumen_conc'
    export_node_field(2, filename + '_conc' + vent + '.exnode', groupname, field_name)
    export_elem_field(filename + '_conc' + vent + '.exelem', groupname, field_name)
    field_name = 'lumen_loss'
    export_node_field(6, filename + '_loss' + vent + '.exnode', groupname, field_name)
    export_elem_field(filename + '_loss' + vent + '.exelem', groupname, field_name)
    field_name = 'sed_loss'
    export_node_field(7, filename + '_sed' + vent + '.exnode', groupname, field_name)
    export_elem_field(filename + '_sed' + vent + '.exelem', groupname, field_name)
    field_name = 'imp_loss'
    export_node_field(8, filename + '_imp' + vent + '.exnode', groupname, field_name)
    export_elem_field(filename + '_imp' + vent + '.exelem', groupname, field_name)
    field_name = 'dif_loss'
    export_node_field(9, filename + '_dif' + vent + '.exnode', groupname, field_name)
    export_elem_field(filename + '_dif' + vent + '.exelem', groupname, field_name)

    field_name = 'bronchial_deposition'
    ne_depos = 13
    export_1d_elem_field(ne_depos, filename + '_depos_field' + vent + '.exelem', groupname, field_name)

    # Export terminal solution
    groupname = 'terminal'
    export_terminal_solution(filename + '_terminal' + vent + '.exnode', groupname)

    # Set the working directory back to it's original location.
    os.chdir(cur_dir)

    # Move the files to the desired directory
    destination_directory = os.path.join(export_directory, str(particle_size))
    shutil.move("airway_result.txt", os.path.join(destination_directory, "airway_result.txt"))
    shutil.move("terminal.txt", os.path.join(destination_directory, "terminal.txt"))


if __name__ == '__main__':
    main()

