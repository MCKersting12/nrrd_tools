import numpy as np
import nrrd
import mcubes
import glob


def main():
    # scale factor moving from bin 2 to microns
    scale_factor = 0.011

    all_files = []
    spacings = []
    mins = []
    lens = []

    # store all files and header info into lists
    i = 0
    for input_file in glob.glob("*.nrrd"):
        f, header = nrrd.read(input_file)
        all_files.append(f)

        spacings.append(header['spacings'])
        mins.append(header['axis mins'])
        lens.append(header['sizes'])
        i += 1

    # find absolute minimum and total x,y,z length needed
    abs_x_min = 9999999999
    abs_y_min = 9999999999
    abs_z_min = 9999999999
    abs_x_max = 0
    abs_y_max = 0
    abs_z_max = 0
    min_x_spacing = 99999999999
    min_y_spacing = 99999999999
    min_z_spacing = 99999999999

    for j in range(len(spacings[:][1])-1):
        print(all_files[j].shape)
        if mins[j][0] < abs_x_min:
            abs_x_min = mins[j][0]

        if mins[j][1] < abs_y_min:
            abs_y_min = mins[j][1]

        if mins[j][2] < abs_z_min:
            abs_z_min = mins[j][2]

        if (mins[j][0]+(lens[j][0]*spacings[j][0])) > abs_x_max:
            abs_x_max = mins[j][0]+(lens[j][0]*spacings[j][0])

        if (mins[j][1]+(lens[j][1]*spacings[j][1])) > abs_y_max:
            abs_y_max = mins[j][1]+(lens[j][1]*spacings[j][1])

        if (mins[j][2]+(lens[j][2]*spacings[j][2])) > abs_z_max:
            abs_z_max = mins[j][2]+(lens[j][2]*spacings[j][2])

        if spacings[j][0] < min_x_spacing:
            min_x_spacing = spacings[j][0]

        if spacings[j][1] < min_y_spacing:
            min_y_spacing = spacings[j][1]

        if spacings[j][2] < min_z_spacing:
            min_z_spacing = spacings[j][2]

    x = int((abs_x_max-abs_x_min)/min_x_spacing)
    y = int((abs_y_max-abs_y_min)/min_y_spacing)
    z = int((abs_z_max-abs_z_min)/min_z_spacing)

    print(x)
    print(y)
    print(z)

    merged_file = np.zeros((x, y, z))

    k = 0
    for file in all_files:

        delta_x = int((abs_x_min - mins[k][0]) / spacings[k][0])
        delta_y = int((abs_y_min - mins[k][1]) / spacings[k][1])
        delta_z = int((abs_z_min - mins[k][2]) / spacings[k][2])

        print('merging file ' + str(k))

        for x in range(lens[k][0]-1):
            if 0 <= x + delta_x < abs_x_max:

                for y in range(lens[k][1]-1):
                    if 0 <= y + delta_y < abs_y_max:

                        for z in range(lens[k][2]-1):
                            if 0 <= z + delta_z < abs_z_max:
                                if file[x+delta_x][y+delta_y][z+delta_z] == 1:

                                    merged_file[x][y][z] = 1

        k += 1

    print("marching cubes...")
    # create mesh
    verts_1, faces_1 = mcubes.marching_cubes(merged_file, 0)

    # scale mesh
    for vert in verts_1:
        vert[0] = ((vert[0]) * spacings[1][1] + abs_x_min) * scale_factor
        vert[1] = ((vert[1]) * spacings[1][2] + abs_y_min) * scale_factor
        vert[2] = ((vert[2]) * spacings[1][3] + abs_z_min) * scale_factor

    mcubes.export_mesh(verts_1, faces_1, 'merged.dae')


def boolean_or(data_1, data_2):
    return (np.logical_or(data_1, data_2))


if __name__ == "__main__":

    main()
