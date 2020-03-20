import numpy as np
import math
import nrrd
from scipy import ndimage
from scipy.ndimage import zoom
import mcubes
import glob
import sys
import os
import code


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
    print("cell body: " + str(header_1['spacings']))
    print("terminal: " + str(header_2['spacings']))

    if x_spacing_1 != 1:
        scale_factor = x_spacing_1
        print("Resampling cell body by factor of " + str(scale_factor))
        cell_body = change_bin(scale_factor, cell_body)
        x_spacing_1, y_spacing_1, z_spacing_1 = 1.0, 1.0, 5.4545
        x_len_1, y_len_1, z_len_1 = x_len_1*scale_factor, y_len_1*scale_factor, z_len_1*scale_factor
    if x_spacing_2 != 1:
        scale_factor = x_spacing_2
        print("Resampling terminal by factor of " + str(scale_factor))
        terminal = change_bin(scale_factor, terminal, z_min=_int(z_min_1), z_max=int(z_min_1+z_len_1*z_spacing_1))
        x_spacing_2, y_spacing_2, z_spacing_2 = 1.0, 1.0, 5.4545
        x_len_2, y_len_2 = x_len_2*scale_factor, y_len_2*scale_factor
        z_len_2 = z_len_1

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

        """
        if delta_x > 0 and delta_y > 0 and delta_z > 0:
            new_terminal[0:x_len_2-delta_x][0:y_len_2-delta_y][0:y_len_2-delta_y] = 
            terminal[delta_x:x_len_2][delta_y:y_len_2][delta_z:z_len_2]
        """

    else:
        new_terminal = terminal

    head, tail = os.path.splitext(terminal_nrrd_file)
    out_file = head + "_ASA.dae"

    print("performing boolean operators...")
    # Remove any part of the cell body that overlaps with the terminal
    # This is to correct for imperfect segmentation, the terminals are typically
    # traced more accurately
    bool_remove_cell = boolean_remove(cell_body, new_terminal)

    dilated_cell = bool_remove_cell

    # dilate more heavily in x and y due to anisotropic voxels
    for z in range(z_len_1):
        dilated_cell[:][:][z] = ndimage.binary_dilation(dilated_cell[:][:][z], iterations=3).astype(terminal.dtype)

    # dilate the cell body to create overlap
    dilated_cell = ndimage.binary_dilation(dilated_cell, iterations=3).astype(int)

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


def change_bin(scale_factor, data, z_min = 0, z_max=2200):

    data = data.astype(bool)
    slices = []
    #Downsample in x and y by given scale_factor
    if z_min == 0 and z_max == 2200:
        for z_slice in range(data.shape[2]):
            current_slice = data[:, :, z_slice]
            zoomed = zoom(current_slice, (scale_factor, scale_factor))
            slices.append(zoomed)
    else:
        for z_slice in range(z_min, z_max):
            current_slice = data[:, :, z_slice]
            zoomed = zoom(current_slice, (scale_factor, scale_factor))
            slices.append(zoomed)

    # cast list of downsampled slices in np array
    new_data = np.asarray(slices)
    #reshuffle array to get back to original order of axes
    new_data = np.swapaxes(new_data, 0, 2)
    new_data = np.swapaxes(new_data, 0, 1)

    return(new_data)

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

    if len(sys.argv) == 3:
        cell_body = sys.argv[1]

        if os.path.isdir(sys.argv[2]):
            for in_file in glob.glob(os.path.join(sys.argv[2], "*.nrrd")):
                if in_file == cell_body:
                    continue
                print(in_file)
                main(cell_body, in_file)
        else:
            main(cell_body, sys.argv[2])
    else:
        print('auto_ASA takes in two arguments: 1- the cell body nrrd file 2- a file or folder containing the inputs')

