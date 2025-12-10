#### File Conversion classes
import os
import numpy as np

# Class for converting to ipfiel file from 1D exelem radius file
class ip_rad:
    def __init__(self,input_file):
        f = open(input_file, "r")
        self.lines = f.readlines()
        f.close()
        self.out_name = input_file.split('.')[0] + '.ipfiel'
        group_name = self.lines[0].split(':')[1]
        self.writing = " CMISS Version 1.21 ipfiel File Version 3\n\
 Heading: "+group_name+"\n"
        num_elem = (self.lines[-3].split(':')[1]).split()[0]
        self.writing += " The number of elements is [  {}]: {}\n\n".format(num_elem,num_elem)
        self.ip_rad = " Element number:        {}\n\
 The field variable value is [ 0.00000D+00]: {:.4g}\n\n"
        
    def write_to_file(self):
        for i in range(8,len(self.lines)-3,3):    
            elem_num = (self.lines[i].split(':')[1]).split()[0]
            elem_val = self.lines[i+2].split()[1]
            self.writing += self.ip_rad.format(int(elem_num), float(elem_val))
        
    def __del__(self):
        w = open(self.out_name, "w")
        w.write(self.writing)
        w.close()
        
# Class for converting to ipelem file from 1D exelem tree
### Usage:
# conversion = ip_elem(input_file.exelem)
# conversion.write_to_file()
# del conversion
class ip_elem:
    def __init__(self,input_file):
        f = open(input_file, "r")
        self.lines = f.readlines()
        f.close()
        self.out_name = input_file.split('.')[0] + '.ipelem'
        group_name = self.lines[0].split(':')[1]
        num_elem = (self.lines[-5].split(':')[1]).split()[0]
        self.writing = " CMISS Version 2.1  ipelem File Version 2\n Heading: {}\n \n".format(group_name) 
        self.writing += " The number of elements is [  {}]: {}\n\n".format(num_elem,num_elem)
        self.ip_elem = " Element number [    {}]:     {}\n\
 The number of geometric Xj-coordinates is [3]: 3\n\
 The basis function type for geometric variable 1 is [1]:  1\n\
 The basis function type for geometric variable 2 is [1]:  1\n\
 The basis function type for geometric variable 3 is [1]:  1\n\
 Enter the {} global numbers for basis 1:"
        
    def ipwriter(self,elem_no,nodes):
        node_no = len(nodes)
        s = self.ip_elem.format(elem_no,elem_no,node_no)
        for i in range(node_no):
            s += '     {}'.format(nodes[i])
        self.writing += s+'\n\n'
        
    def write_to_file(self):        
        for i in range(31,len(self.lines)-3,5):    
            elem_num = (self.lines[i].split(':')[1]).split()[0]
            elem_val = self.lines[i+2].split()
            self.ipwriter(int(elem_num),elem_val)
        
    def __del__(self):
        w = open(self.out_name, "w")
        w.write(self.writing)
        w.close()


# Class for converting to ipnode file from 1D exnode tree
### Usage:
# conversion = ip_node(input_file.exnode)
# conversion.write_to_file()
# del conversion
# class ip_node:
#     def __init__(self,input_file):
#         f = open(input_file, "r")
#         self.lines = f.readlines()
#         f.close()
#         self.out_name = input_file.split('.')[0] + '.ipnode' # should be ipnode, not ipfiel
#         group_name = self.lines[0].split(':')[1]
#
#         for index, line in enumerate(self.lines):
#             if 'Node:' in line:
#                 self.start = index
#                 second_line = index+2
#                 break
#         if 'Node:' in self.lines[second_line]:
#             self.case = 1
#             self.end = len(self.lines)-2
#             self.skip = second_line - self.start
#             num_node = (self.lines[self.end].split(':')[1]).strip()
#
#         else:
#             # This Variant needs to be double checked
#             self.case = 0
#             self.end = len(self.lines)-4
#             self.skip = 4
#             num_node = (self.lines[self.end].split(':')[1]).strip()
#
#         self.writing = " CMISS Version 2.1  ipnode File Version 2\n Heading: {}\n ".format(group_name)
#         self.writing += " The number of nodes is [  {}]: {}\n".format(num_node,num_node)
#         self.writing += " Number of coordinates [ 3]:  3\n\
#  Do you want prompting for different versions of nj=1 [N]? n\n\
#  Do you want prompting for different versions of nj=2 [N]? n\n\
#  Do you want prompting for different versions of nj=3 [N]? n\n\
#  The number of derivatives for coordinate 1 is [0]:  \n\
#  The number of derivatives for coordinate 2 is [0]:  \n\
#  The number of derivatives for coordinate 3 is [0]:  \n\n"
#
#     def ipwriter(self,node_curr,node_val):
#         self.writing += ' Node number [ {}]:  {}\n\
#  The Xj(1) coordinate is [ 0.00000E+00]:  {} \n\
#  The Xj(2) coordinate is [ 0.00000E+00]:  {} \n\
#  The Xj(3) coordinate is [ 0.00000E+00]:  {} \n\n'.format(node_curr,node_curr, node_val[0], node_val[1], node_val[2])
#
#     def write_to_file(self):
#         for i in range(self.start,self.end+self.skip,self.skip):    # This end may not work with case ==0
#             node_num = self.lines[i].split(':')[1]
#             if self.case:
#                 temp= self.lines[i+1].split()
#                 coord_vals = [float(x) for x in temp]
#             else:
#                 coord_vals = [float(self.lines[i+j]) for j in range(1,4)]
#             self.ipwriter(int(node_num),coord_vals)
#
#     def __del__(self):
#         w = open(self.out_name, "w")
#         w.write(self.writing)
#         w.close()

class ip_node:
    def __init__(self, input_file):
        with open(input_file, "r") as f:
            self.lines = f.readlines()

        self.out_name = input_file.split('.')[0] + '.ipnode'

        # Robustly find the number of nodes
        num_node = None
        for line in reversed(self.lines):
            if line.strip().startswith('Node:'):
                num_node = int(line.split()[-1])
                break
        if num_node is None:
            raise ValueError("Could not find number of nodes in the input file.")

        group_name = "Airway_Node_Geometry"  # or appropriate name
        self.writing = f" CMISS Version 2.1  ipnode File Version 2\n Heading: {group_name}\n"
        self.writing += f" The number of nodes is [ {num_node}]: {num_node}\n"
        self.writing += (" Number of coordinates [ 3]:  3\n"
                         " Do you want prompting for different versions of nj=1 [N]? n\n"
                         " The number of derivatives for geometric variables [0]: 0\n\n")

    def ipwriter(self, node_curr, node_val):
        self.writing += f' Node number [ {node_curr}]:  {node_curr}\n'
        self.writing += f' The Xj(1) coordinate is [ 0.00000E+00]: {node_val[0]}\n'
        self.writing += f' The Xj(2) coordinate is [ 0.00000E+00]: {node_val[1]}\n'
        self.writing += f' The Xj(3) coordinate is [ 0.00000E+00]: {node_val[2]}\n\n'

    def write_to_file(self):
        node_num = None
        i = 0
        while i < len(self.lines):
            line = self.lines[i]
            if 'Node:' in line:
                node_num = int(line.strip().split()[1])
                x = float(self.lines[i + 1].strip())
                y = float(self.lines[i + 2].strip())
                z = float(self.lines[i + 3].strip())
                self.ipwriter(node_num, [x, y, z])
                i += 4
            else:
                i += 1

        with open(self.out_name, "w") as w:
            w.write(self.writing)


if __name__ == "__main__":
    # f1 = "geometries/grown-007.exnode"
    # f2 = "geometries/grown-007.exelem"
    # f3 = "geometries/grown-007_radius.exelem"
    #
    f1 = "../../Packages/lung-group-examples/growing_tri_surface_Tawhai2023/output/004_Airway_Full.exnode"
    f2 = "../../Packages/lung-group-examples/growing_tri_surface_Tawhai2023/output/004_Airway_Full.exelem"
    #f3 = "geometries/grown-007_radius.exelem"
    conversion = ip_node(f1)
    conversion.write_to_file()
    del conversion
    conversion = ip_elem(f2)
    conversion.write_to_file()
    del conversion
    # conversion = ip_rad(f3)
    # conversion.write_to_file()
    # del conversion



