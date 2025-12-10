#!/usr/bin/env python

import os
import argparse
from aether.diagnostics import set_diagnostics_on
from aether.indices import define_problem_type, ventilation_indices, get_ne_radius
from aether.geometry import define_node_geometry, define_1d_elements, define_rad_from_file, define_rad_from_geom
from aether.geometry import append_units, refine_1d_elements, reorder_tree
from aether.geometry import list_tree, scale_tree, scale_radii
from aether.ventilation import evaluate_vent
from aether.exports import export_1d_elem_geometry, export_node_geometry, export_elem_field, export_1d_elem_field, \
    export_terminal_solution


# -----------------------------------------------------------------------------------------------------


def HLA_subject_data():
    # dictionary to contain the subject lung measurement data
    # TJ - for H1335, Rtrach = hght * 0.74394657607831 + 1.63* 3.1634738580296 = 6.44348996520372
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

#
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
    parser = argparse.ArgumentParser(description='Run ventilation for a selected subject')
    parser.add_argument('-s', '--subject', help='Input the subject name', required=True)

    args = parser.parse_args()
    subject_data = HLA_subject_data()
    abbrv_subject = args.subject
    for key in subject_data:
        if abbrv_subject in key:
            subject = key

    volume = 'FRC'

    set_diagnostics_on(False)

    # Read settings
    define_problem_type('ventilation')
    ventilation_indices()

    #export_directory = 'output'
    export_directory = os.path.join('output', str(abbrv_subject))

    if not os.path.exists(export_directory):
        os.makedirs(export_directory)

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
    define_rad_from_file('geometries/upper_template.ipfiel', 1.0) # param to scale deadspace
    scale_trachea = trachea_rad / 8.0  # divide by the trachea size in the template file
    scale_radii(scale_trachea)  # scale all of the template by the same value

    min_rad = 0.25 # adjust to make L/D between [2.8, 3.0]
    define_rad_from_geom('fit', 0.0, 'inlet', min_rad)
    list_tree()

    refine_elements = list(range(1,21))
    refine_1d_elements(refine_elements, 10)
    reorder_tree()

    append_units()

    # Set the working directory to this file's directory and then reset after running simulation.
    file_location = os.path.dirname(os.path.abspath(__file__))
    cur_dir = os.getcwd()
    os.chdir(file_location)

    # Run simulation.
    evaluate_vent()

    # Set the working directory back to it's original location.
    os.chdir(cur_dir)

    # Output results

    # Export airway nodes and elements
    group_name = 'vent_model'
    export_filename = export_directory + '/' + subject + '_geometry'
    export_1d_elem_geometry(export_filename, group_name)
    export_node_geometry(export_filename, group_name)
    field_name = 'flow'
    export_filename = export_directory + '/' + subject + '_ventilation_field'
    export_1d_elem_field(6, export_filename, group_name, field_name)

    # Export element field for radius
    ne_radius = get_ne_radius()
    field_name = 'radius'
    export_filename = export_directory + '/' + subject + '_radius_field'
    export_1d_elem_field(ne_radius, export_filename, group_name, field_name)

    # Export element field for resistance
    field_name = 'resist'
    export_filename = export_directory + '/' + subject + '_resist_field'
    export_1d_elem_field(4, export_filename, group_name, field_name)

    # Export element field for total resistance
    field_name = 't_resist'
    export_filename = export_directory + '/' + subject + '_t_resist_field'
    export_1d_elem_field(5, export_filename, group_name, field_name)

    # Export element field for pressure drop
    field_name = 'press_drop'
    export_filename = export_directory + '/' + subject + '_pressure_field'
    export_1d_elem_field(2, export_filename, group_name, field_name)

    # Export element field for branched ventilation field
    field_name = 'vent_avd'
    export_filename = export_directory + '/' + subject + '_ventavd_field'
    export_1d_elem_field(7, export_filename, group_name, field_name)

    # Export terminal solution
    group_name = 'terminal'
    export_filename = export_directory + '/' + subject + '_terminal.exnode'
    export_terminal_solution(export_filename, group_name)


if __name__ == '__main__':
    main()