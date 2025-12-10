import os
import argparse
import subprocess
import re
import math
from aether.diagnostics import set_diagnostics_on
from aether.indices import define_problem_type, ventilation_indices, get_ne_radius
from aether.geometry import define_node_geometry, define_1d_elements, define_rad_from_geom, append_units, list_tree, \
    evaluate_ordering, group_elem_by_parent
from aether.ventilation import evaluate_vent
from aether.exports import export_1d_elem_geometry, export_node_geometry, export_1d_elem_field, export_terminal_solution


#     '004': {'h_ratio': 1.165, 'trachea_rad': 6.1225, 'min_rad': 0.13, 'FRC': 1.14},
#     '005': {'h_ratio': 1.145, 'trachea_rad': 6.70875, 'min_rad': 0.13, 'FRC': 1.632},
#     '007': {'h_ratio': 1.15, 'trachea_rad': 6.3625, 'min_rad': 0.12, 'FRC': 1.434},
#     '009': {'h_ratio': 1.14, 'trachea_rad': 6.150625, 'min_rad': 0.15, 'FRC': 1.622},
#     '012': {'h_ratio': 1.15, 'trachea_rad': 6.155, 'min_rad': 0.15, 'FRC': 1.402},
#     '013': {'h_ratio': 1.17, 'trachea_rad': 8.868, 'min_rad': 0.16, 'FRC': 1.651},
#     '021': {'h_ratio': 1.16, 'trachea_rad': 6.486, 'min_rad': 0.16, 'FRC': 1.302},
#     '025': {'h_ratio': 1.1475, 'trachea_rad': 6.765, 'min_rad': 0.13, 'FRC': 1.407},


# Dictionary containing subject-specific values
SUBJECT_PARAMETERS = {
    '004': {'h_ratio': 1.1525, 'trachea_rad': 4.898, 'min_rad': 0.13, 'FRC': 1.14},
        '005': {'h_ratio': 1.1275, 'trachea_rad': 4.4725, 'min_rad': 0.13, 'FRC': 1.632},
        '007': {'h_ratio': 1.14, 'trachea_rad': 5.09, 'min_rad': 0.12, 'FRC': 1.434},
        '009': {'h_ratio': 1.13, 'trachea_rad': 4.9205, 'min_rad': 0.15, 'FRC': 1.622},
        '012': {'h_ratio': 1.15, 'trachea_rad': 6.155, 'min_rad': 0.15, 'FRC': 1.402},
        '013': {'h_ratio': 1.175, 'trachea_rad': 8.868, 'min_rad': 0.16, 'FRC': 1.651},
        '021': {'h_ratio': 1.15, 'trachea_rad': 5.405, 'min_rad': 0.16, 'FRC': 1.302},
        '025': {'h_ratio': 1.13, 'trachea_rad': 4.51, 'min_rad': 0.13, 'FRC': 1.407},
    '026': {'h_ratio': 1.1125, 'trachea_rad': 3.0535, 'min_rad': 0.13, 'FRC': 1.266},
    '029': {'h_ratio': 1.1475, 'trachea_rad': 4.748, 'min_rad': 0.15, 'FRC': 1.681}
}


def write_parameters(subject, FRC):
    """Writes FRC and other ventilation parameters to 'params_evaluate_flow.txt'."""
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
        opf.write("volume_target 500000\n")
        opf.write("pmus_step  -196.133\n")
        opf.write("expiration_type active \n")
        opf.write("chest_wall_compliance 2039.4324\n")

    print(f"Parameters file updated: {params_path}")




def main():
    parser = argparse.ArgumentParser(description='Run ventilation for a selected subject')
    parser.add_argument('-s', '--subject', help='Input the subject name', required=True)

    args = parser.parse_args()
    subject = args.subject

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

    # Write parameters file with correct FRC
    write_parameters(subject, FRC)

    set_diagnostics_on(False)

    # Read settings
    define_problem_type('ventilation')
    ventilation_indices()

    filename = f'{subject}-virtual-500'
    export_directory = os.path.join('output', filename)
    os.makedirs(export_directory, exist_ok=True)


    # Define geometry based on subject
    airway_geometry_path = f'../../Packages/lung-group-examples/growing_tri_surface_Tawhai2023/output/{subject}_Airway_virtual'
    define_node_geometry(airway_geometry_path)
    define_1d_elements(airway_geometry_path)

    define_rad_from_geom('horsf', h_ratio, 'inlet', trachea_rad, 'all')

    list_tree()
    append_units()

    # Run ventilation simulation
    evaluate_vent()

    # Export results
    group_name = 'vent_model'

    export_filename = os.path.join(export_directory, f"{subject}_geometry")
    export_1d_elem_geometry(export_filename, group_name)
    export_node_geometry(export_filename, group_name)

    # Export ventilation fields
    field_name = 'flow'
    export_filename = os.path.join(export_directory, f"{subject}_ventilation_field")
    export_1d_elem_field(6, export_filename, group_name, field_name)

    # Export element fields
    field_mappings = {
        'radius': get_ne_radius(),
        'resist': 4,
        't_resist': 5,
        'press_drop': 2,
        'vent_avd': 7,
        'surface_area': 11
    }

    for field_name, index in field_mappings.items():
        export_filename = os.path.join(export_directory, f"{subject}_{field_name}_field")
        export_1d_elem_field(index, export_filename, group_name, field_name)

    # Export terminal solution
    export_filename = os.path.join(export_directory, f"{subject}_terminal.exnode")
    export_terminal_solution(export_filename, 'terminal')

    print(f"Ventilation simulation completed for subject {subject}.")


if __name__ == '__main__':
    main()
