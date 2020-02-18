import nrrd
import kimimaro
from scipy.ndimage import zoom
import numpy as np
import os
import pySWC
import glob
import sys


def main(input_file):

    labels, header = nrrd.read(input_file)
    print("Original shape: " + str(labels.shape))
    mins = header['axis mins']

    labels = labels.astype(bool)
    slices = []
    #Downsample in x and y by given scale_factor
    scale_factor = 0.5
    for z_slice in range(labels.shape[2]):
        current_slice = labels[:, :, z_slice]
        zoomed = zoom(current_slice, (scale_factor, scale_factor))
        slices.append(zoomed)

    # cast list of downsampled slices in np array
    labels = np.asarray(slices)
    #reshuffle array to get back to original order of axes
    labels = np.swapaxes(labels, 0, 2)
    labels = np.swapaxes(labels, 0, 1)

    skels = kimimaro.skeletonize(labels,
        teasar_params={
            'scale': 12,
            'const': 0, # physical units
            'pdrf_exponent': 4,
            'pdrf_scale': 100000,
            'soma_detection_threshold': 1100, # physical units
            'soma_acceptance_threshold': 3500, # physical units
            'soma_invalidation_scale': 1.0,
            'soma_invalidation_const': 300, # physical units
            'max_paths': 50, # default None
        },
        # object_ids=[ ... ], # process only the specified labels
        # extra_targets_before=[ (27,33,100), (44,45,46) ], # target points in voxels
        # extra_targets_after=[ (27,33,100), (44,45,46) ], # target points in voxels
        dust_threshold=1000, # skip connected components with fewer than this many voxels
        anisotropy=(1,1,1), # default True
        fix_branching=True, # default True
        fix_borders=True, # default True
        progress=True, # default False, show progress bar
        parallel=1, # <= 0 all cpu, 1 single process, 2+ multiprocess
        parallel_chunk_size=100, # how many skeletons to process before updating progress bar
    )

    name, ext = os.path.splitext(input_file)

    with open(name + ".swc", "wt") as f:
        f.write(skels[1].to_swc())
        swc = pySWC.Swc(name + ".swc")

        #Translate to nrrd minimums
        swc.translate(mins[0], mins[1], mins[2])

        #correct for scaling of x and y
        swc.anisotropic_scale(1/scale_factor, 1/scale_factor, 1)

        #Move into micron units
        swc.scale(0.011)
        swc.save_file(name + "_scaled.swc")


if __name__ == "__main__":

    if len(sys.argv) > 1:

        if os.path.isdir(sys.argv[1]):
            for in_file in glob.glob(os.path.join(sys.argv[1], "*.nrrd")):
                main(in_file)
        else:
            main(sys.argv[1])

    else:
        print("Usage: python downsample_teaser.py [folder containing nrrds]")