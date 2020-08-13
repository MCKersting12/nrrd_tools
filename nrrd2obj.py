"""
Converts a volume image file (in nrrd format) to a mesh file (in obj format)
"""

import nrrd
import mcubes
import sys
import os


def main(input_file):

    print('loading file...')
    input_nrrd, header = nrrd.read(input_file)

    x_spacing, y_spacing, z_spacing = header['spacings']
    x_min, y_min, z_min = header['axis mins']
    x_len, y_len, z_len = header['sizes']

    head, tail = os.path.splitext(input_file)
    out_file = head + ".obj"

    print("marching cubes...")
    # create mesh
    verts, faces = mcubes.marching_cubes(input_nrrd, 0)
    
    for vert in verts:
        vert[0] = ((vert[0])* x_spacing + x_min)
        vert[1] = ((vert[1])* y_spacing + y_min)
        vert[2] = ((vert[2])* z_spacing + z_min)

    mcubes.export_obj(verts, faces, out_file)


if __name__ == "__main__":

    if len(sys.argv) > 1:
        main(sys.argv[1])

    else:
        print('nrrd2dae.py takes in one argument: the nrrd file to be converted')