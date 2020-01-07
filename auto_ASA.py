import numpy as np
import math
import nrrd
from scipy import ndimage
import mcubes
import glob
import sys
import os


def main(cell_nrrd_file, terminal_nrrd_file):

    # scale factor moving from bin 2 to microns
    scale_factor = 0.011

    # grab both nrrd files and read them into numpy arrays
    cell_body, header_1 = nrrd.read(cell_nrrd_file)
    terminal, header_2 = nrrd.read(terminal_nrrd_file)

    x_spacing_1, y_spacing_1, z_spacing_1 = header_1['spacings']
    x_min_1, y_min_1, z_min_1 = header_1['axis mins']
    x_len_1, y_len_1, z_len_1 = header_1['sizes']
    x_spacing_2, y_spacing_2, z_spacing_2 = header_2['spacings']
    x_min_2, y_min_2, z_min_2 = header_2['axis mins']
    x_len_2, y_len_2, z_len_2 = header_2['sizes']

    # if the dimensions are not the same then we resample the terminal
    if x_spacing_1 != x_spacing_2 or x_min_1 != x_min_2 or y_min_1 != y_min_2 \
        or y_spacing_1 != y_spacing_2 or z_min_1 != z_min_2 or \
            z_spacing_1 != z_spacing_2:

        delta_x = int((x_min_1 - x_min_2) / x_spacing_1)
        delta_y = int((y_min_1 - y_min_2) / y_spacing_1)
        delta_z = int((z_min_1 - z_min_2) / z_spacing_1)

        print("resampling terminal file...")
        new_terminal = np.zeros((x_len_1, y_len_1, z_len_1))

        for x in range(x_len_1):
            if 0 <= x+delta_x < x_len_2:

                for y in range(y_len_1):
                    if 0 <= y + delta_y < y_len_2:

                        for z in range(z_len_1):
                            if 0 <= z + delta_z < z_len_2:

                                new_terminal[x][y][z] = terminal[x+delta_x][y+delta_y][z+delta_z]

        # should be able to resample the terminal based in indexing, looping is inefficient
        # new_terminal[:int(x_max_2-x_min_1)][:int(y_max_2-y_min_1)][:int(z_max_2-z_min_2)] =
        # terminal[delta_x:][delta_y:][delta_z:]

    else:
        new_terminal = terminal

    head, tail = os.path.splitext(terminal_nrrd_file)
    out_file = head + "_ASA.dae"

    print("performing boolean operators...")
    # Remove any part of the cell body that overlaps with the terminal
    # This is to correct for imperfect segmentations, the terminals are typically
    # traced more accurately
    bool_remove_cell = boolean_remove(cell_body, new_terminal)

    dilated_cell = bool_remove_cell

    # dilate more heavily in x and y due to anisotropic voxels
    for z in range(z_len_1):
        dilated_cell[:][:][z] = ndimage.binary_dilation(dilated_cell[:][:][z], iterations=3).astype(terminal.dtype)

    # dilate the cell body to create overlap
    dilated_cell = ndimage.binary_dilation(bool_remove_cell, iterations=3).astype(int)

    # grab the overlapping area
    asa_matrix = boolean_and(dilated_cell, new_terminal)

    print("marching cubes...")
    # create mesh
    verts_1, faces_1 = mcubes.marching_cubes(asa_matrix, 0)

    # scale mesh
    for vert in verts_1:
        vert[0] = ((vert[0])*x_spacing_1 + x_min_1)*scale_factor
        vert[1] = ((vert[1])*y_spacing_1 + y_min_1)*scale_factor
        vert[2] = ((vert[2])*z_spacing_1 + z_min_1)*scale_factor

    mcubes.export_mesh(verts_1, faces_1, out_file)
    surface_area = calculate_surface_area(verts_1, faces_1)
    print(terminal_nrrd_file + " apposed surface area = " + str(surface_area/2.0))
    print('...')


def boolean_remove(data_1, data_2):

    return(np.greater(data_1, data_2))


def boolean_and(data_1, data_2):

    return(np.logical_and(data_1, data_2))


def boolean_or(data_1, data_2):

    return(np.logical_or(data_1, data_2))


def calculate_surface_area(verts, faces):
    surface_area = 0.0
    print("finding surface area...")
    for face in faces:
        vert_1, vert_2, vert_3 = face

        vert_1x, vert_1y, vert_1z = verts[vert_1]
        vert_2x, vert_2y, vert_2z = verts[vert_2]
        vert_3x, vert_3y, vert_3z = verts[vert_3]

        dis1_2_x = vert_1x - vert_2x
        dis1_2_y = vert_1y - vert_2y
        dis1_2_z = vert_1z - vert_2z
        dis1_3_x = vert_1x - vert_3x
        dis1_3_y = vert_1y - vert_3y
        dis1_3_z = vert_1z - vert_3z
        
        dis2_1_x = vert_2x - vert_1x
        dis2_1_y = vert_2y - vert_1y
        dis2_1_z = vert_2z - vert_1z
        dis3_1_x = vert_3x - vert_1x
        dis3_1_y = vert_3y - vert_1y
        dis3_1_z = vert_3z - vert_1z

        # Carol starts here
        b = math.sqrt((dis1_2_x**2)+(dis1_2_y**2)+(dis1_2_z**2))
        d = math.sqrt((dis1_3_x**2)+(dis1_3_y**2)+(dis1_3_z**2))

        if b*d != 0.0:
            theta = math.acos(round(((dis3_1_x*dis2_1_x + dis3_1_y*dis2_1_y + dis3_1_z*dis2_1_z) / (b*d)), 4))

            # check for nans
            if math.isnan(theta):
                theta = 0.0

            area = b*d*math.sin(theta) / 2.0

            if math.isnan(area):
                area = 0.0

            surface_area += area

    return(surface_area)


if __name__ == "__main__":

    if len(sys.argv) > 1:
        cell_file = sys.argv[1]

        for input_file in glob.glob("*.nrrd"):
            print(input_file)

            if input_file != cell_file:
                main(cell_file, input_file)

    else:
        print('auto_ASA.py takes in one argument: the cell body name (or just the object being contacted)')
        print('all files need to be in the same directory as the script and cell file')
