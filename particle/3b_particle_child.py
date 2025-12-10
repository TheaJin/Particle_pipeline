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
    define_rad_from_geom, refine_1d_elements, reorder_tree, set_initial_volume, apply_cluster_constriction
from aether.geometry import list_tree, scale_tree, scale_radii
from aether.particle_transport import solve_particles_decoupled
from aether.ventilation import evaluate_uniform_flow
from aether.field_utilities import scale_flow_to_inlet

# Dictionary containing subject-specific values
# SUBJECT_PARAMETERS = {
#     '004': {'h_ratio': 1.1525, 'trachea_rad': 4.898, 'min_rad': 0.13, 'FRC': 1.14},
#         '005': {'h_ratio': 1.1275, 'trachea_rad': 4.4725, 'min_rad': 0.13, 'FRC': 1.632},
#         '007': {'h_ratio': 1.14, 'trachea_rad': 5.09, 'min_rad': 0.12, 'FRC': 1.434},
#         '009': {'h_ratio': 1.13, 'trachea_rad': 4.9205, 'min_rad': 0.15, 'FRC': 1.622},
#         '012': {'h_ratio': 1.15, 'trachea_rad': 6.155, 'min_rad': 0.15, 'FRC': 1.402},
#         '013': {'h_ratio': 1.175, 'trachea_rad': 8.868, 'min_rad': 0.16, 'FRC': 1.651},
#         '021': {'h_ratio': 1.15, 'trachea_rad': 5.405, 'min_rad': 0.16, 'FRC': 1.302},
#         '025': {'h_ratio': 1.13, 'trachea_rad': 4.51, 'min_rad': 0.13, 'FRC': 1.407},
#     '026': {'h_ratio': 1.1125, 'trachea_rad': 3.0535, 'min_rad': 0.13, 'FRC': 1.266},
#     '029': {'h_ratio': 1.1475, 'trachea_rad': 4.748, 'min_rad': 0.15, 'FRC': 1.681}
# }
SUBJECT_PARAMETERS = {
    '004': {'h_ratio': 1.165, 'trachea_rad': 6.1225, 'min_rad': 0.13, 'FRC': 1.14},
    '005': {'h_ratio': 1.145, 'trachea_rad': 6.70875, 'min_rad': 0.13, 'FRC': 1.632},
    '007': {'h_ratio': 1.15, 'trachea_rad': 6.3625, 'min_rad': 0.12, 'FRC': 1.434},
    '009': {'h_ratio': 1.14, 'trachea_rad': 6.150625, 'min_rad': 0.15, 'FRC': 1.622},
    '012': {'h_ratio': 1.15, 'trachea_rad': 6.155, 'min_rad': 0.15, 'FRC': 1.402},
    '013': {'h_ratio': 1.17, 'trachea_rad': 8.868, 'min_rad': 0.16, 'FRC': 1.651},
    '021': {'h_ratio': 1.16, 'trachea_rad': 6.486, 'min_rad': 0.16, 'FRC': 1.302},
    '025': {'h_ratio': 1.1475, 'trachea_rad': 6.765, 'min_rad': 0.13, 'FRC': 1.407},
    '026': {'h_ratio': 1.1125, 'trachea_rad': 3.0535, 'min_rad': 0.13, 'FRC': 1.266},
    '029': {'h_ratio': 1.1475, 'trachea_rad': 4.748, 'min_rad': 0.15, 'FRC': 1.681}
}

def write_parameters(subject, FRC):
    """Writes FRC and other parameters to 'params_evaluate_flow.txt'."""
    params_path = 'Parameters/params_evaluate_flow.txt'

    with open(params_path, "w") as opf:
        opf.write("FRC        %4.2f\n" % FRC)
        opf.write("constrict  1.0\n")
        opf.write("T_interval 4.0\n")
        opf.write("Gdirn      3\n")
        opf.write("press_in   0.0 \n")
        opf.write("COV        0.1\n")

        # oneG upright
        opf.write("RMaxMean   1.29\n")
        opf.write("RMinMean   0.79\n")

        opf.write("i_to_e_ratio 1.0\n")
        opf.write("refvol     0.5\n")
        opf.write("volume_target 200000\n")
        opf.write("pmus_step  -196.133\n")
        opf.write("expiration_type active \n")
        opf.write("chest_wall_compliance 2039.4324\n")

    print(f"Parameters file updated: {params_path}")


def main():
    parser = argparse.ArgumentParser(description='Run particle transport for a selected subject')
    parser.add_argument('-s', '--subject', help='Input the subject name', required=True)
    parser.add_argument('-p', '--particle_sizes', nargs='+', help='Input particle sizes in microns (space-separated)', required=True)
    parser.add_argument('-g', '--gravity', help='Gravity direction (0/1/2/3/-3)', required=True)
    parser.add_argument('-v', '--vent', help='Input ventilation file code', required=False, default='')
    parser.add_argument('-c', '--constriction_ratio', type=float, default=0.8, help='Homogenous constriction param')
    parser.add_argument('--cluster', action='store_true',
                        help='Apply clustered bronchoconstriction using built-in Fortran subroutine')

    args = parser.parse_args()
    subject = args.subject
    particle_sizes = [float(p) for p in args.particle_sizes]
    Gdirn = int(args.gravity)
    vent = args.vent
    volume = 'FRC'
    constriction_ratio = float(args.constriction_ratio)

    # Retrieve subject-specific values
    if subject in SUBJECT_PARAMETERS:
        params = SUBJECT_PARAMETERS[subject]
        h_ratio = params['h_ratio']
        trachea_rad = params['trachea_rad']
        min_rad = params['min_rad']
        FRC = params['FRC']
    else:
        print(f"Error: No parameters found for subject {subject}. Check your subject ID.")
        return

    # Write the parameters file with the correct FRC
    write_parameters(subject, FRC)

    for particle_size in particle_sizes:
        set_diagnostics_on(False)

        export_directory = os.path.join('output', subject + '-constrict')
        print(export_directory)
        if not os.path.exists(export_directory):
            os.makedirs(export_directory)
        plot_me = False

        define_problem_type('particle_transport')

        airway_geometry_path = f'../../Packages/lung-group-examples/growing_tri_surface_Tawhai2023/output/{subject}_Airway_virtual'
        define_node_geometry(airway_geometry_path)
        define_1d_elements(airway_geometry_path)

        define_rad_from_geom('horsf', h_ratio, 'inlet', trachea_rad, 'all')

        # min_allowed_radius = 0.1
        # scale_radii(constriction_ratio)  # scale all of the template by the same value
        # define_rad_from_geom('fit', 0.0, 'inlet', min_allowed_radius)

        # 🧠 Clustered bronchoconstriction (calls your Fortran!)
        if args.cluster:
            # 12 is an example parent element — you can adjust based on subject or order
            apply_cluster_constriction(12)
            print("✔️ Clustered bronchoconstriction applied via Fortran subroutine.")

        else:
            # 🔁 Homogeneous constriction as default
            scale_radii(constriction_ratio)
            print(f"✔️ Homogeneous constriction applied with scale factor {constriction_ratio}")

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
        # TJ - be careful what to read!!!!!!!!!
        import_filename = export_directory + '/' + subject + '_ventilation_field' + vent
        import_ventilation(import_filename + '.exelem')

        # this is a normalised field so needs to be scaled to our required flow value
        tidal_volume = 200000  # mm3
        time_inspiration = 2.0  # seconds
        flow_rate = tidal_volume / time_inspiration

        # scale flow to alveolar ventilation
        scale_flow_to_inlet(flow_rate, 'V')

        initial_concentration = 0
        inlet_concentration = 0.001  # consistent with inner workings of lungsim

        # solve particle transport
        solve_particles_decoupled(initial_concentration, inlet_concentration, particle_size)

        # name = "{0}_p{1}".format(subject, int(particle_size))
        # filename = os.path.join(export_directory, str(particle_size), name)
        name = f"{subject}_p{particle_size:.2f}"
        filename = os.path.join(export_directory, f"{particle_size:.2f}", name)
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
        #destination_directory = os.path.join(export_directory, str(particle_size))

        destination_directory = os.path.join(export_directory, f"{particle_size:.2f}")

        # Explicitly ensure directory exists:
        os.makedirs(destination_directory, exist_ok=True)

        shutil.move("airway_result.txt", os.path.join(destination_directory, "airway_result.txt"))
        shutil.move("terminal.txt", os.path.join(destination_directory, "terminal.txt"))
        shutil.move("final_results.txt", os.path.join(destination_directory, "final_results.txt"))

if __name__ == '__main__':
    main()
