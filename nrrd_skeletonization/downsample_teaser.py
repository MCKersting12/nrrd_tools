import nrrd
import kimimaro
from scipy.ndimage import zoom
import numpy as np
import code
import os
import pySWC
import glob

for input_file in glob.glob("*.nrrd"):
    labels, header = nrrd.read(input_file)
    print("Original shape: " + str(labels.shape))
    mins = header['axis mins']

    labels = labels.astype(bool)
    xx,yy,zz = labels.shape
    slices = []
    scale_factor = 0.5
    for z_slice in range(labels.shape[2]):
        current_slice = labels[:, :, z_slice]
        zoomed = zoom(current_slice, (scale_factor, scale_factor))
        slices.append(zoomed)

    labels = np.asarray(slices)
    labels = np.swapaxes(labels, 0, 2)
    labels = np.swapaxes(labels, 0, 1)
    print("New shape: " + str(labels.shape))
    """
    scales = [1, 2, 4, 8, 16]
    consts = [1, 200, 500, 1000]
    scale = 16 const = 1 was best result
    """
    """
    scales = [12, 16, 20, 24]
    consts = [0, 1, 2]
    Scale = 12 const = 0 was best result
    """

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
        #code.interact(local=locals())
        f.write(skels[1].to_swc())
        swc = pySWC.Swc(name + ".swc")
        #Translate to nrrd minimums
        swc.translate(mins[0], mins[1], mins[2])
        #correct for scaling of x and y
        swc.anisotropic_scale(1/scale_factor, 1/scale_factor, 1)
        #Move into micron units
        swc.scale(0.011)
        swc.save_file(name + "_scaled.swc")
