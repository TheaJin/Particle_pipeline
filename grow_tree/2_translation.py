# import nrrd
import os
import numpy as np


def read_node_file(filename):
    """
    Reads the file and extracts node data.
    Returns a dictionary with node numbers as keys and (x, y, z) coordinates as values.
    """
    nodes = {}
    with open(filename, 'r') as file:
        lines = file.readlines()

        for i in range(len(lines)):
            line = lines[i].strip()

            # Check if we reached a line starting with "Node:"
            if line.startswith("Node:"):
                node_number = int(line.split()[1])
                # Extract x, y, z coordinates from the following three lines
                x = float(lines[i + 1].strip())
                y = float(lines[i + 2].strip())
                z = float(lines[i + 3].strip())
                nodes[node_number] = (x, y, z)

    return nodes


def translate_nodes(nodes, translation_vector):
    """
    Applies the translation vector to each node's coordinates.
    """
    translated_nodes = {}
    tx, ty, tz = translation_vector
    for node_number, (x, y, z) in nodes.items():
        # Apply translation vector to each coordinate
        scaling = 1
        translated_nodes[node_number] = (scaling* x + tx, scaling* y + ty, scaling* z + tz)
    return translated_nodes


def write_translated_nodes(output_filename, translated_nodes):
    """
    Writes the translated nodes to a new file.
    """
    with open(output_filename, 'w') as outfile:
        outfile.write(' Group name : MAC\n')
        outfile.write(' #Fields=1\n')
        outfile.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        outfile.write('   x.  Value index= 1, #Derivatives= 0\n')
        outfile.write('   y.  Value index= 2, #Derivatives= 0\n')
        outfile.write('   z.  Value index= 3, #Derivatives= 0\n')
        for node_number, (x, y, z) in translated_nodes.items():
            outfile.write(f" Node: {node_number}\n")
            outfile.write(f" {x: .15e}\n")
            outfile.write(f" {y: .15e}\n")
            outfile.write(f" {z: .15e}\n")


# Read the NRRD segmentation file
root_folder = '/hpc/gjin250/'

# origin = np.array([-5.00, -32.00, 250.00]) # 004
origin = np.array([5.00, -5.00, -17.00])

# Main function to run the translation
def translate_node_file(input_filename, output_filename, translation_vector):
    nodes = read_node_file(input_filename)
    translated_nodes = translate_nodes(nodes, translation_vector)
    write_translated_nodes(output_filename, translated_nodes)


file_path = os.path.join(root_folder, 'Packages', 'lung-group-examples', 'geometry', 'template_009.exnode')
# Set parameters for translation
output_filename = os.path.join(root_folder, 'Packages', 'lung-group-examples', 'geometry', 'template_026_translated.exnode')

# Run the translation and write output
translate_node_file(file_path, output_filename, origin)