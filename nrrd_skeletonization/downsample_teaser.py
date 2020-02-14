import nrrd
import kimimaro
from scipy.ndimage import zoom
import numpy as np
import code

labels, header = nrrd.read('VCN_c06_dendrite01_part1.nrrd')
print("Original shape: " + str(labels.shape))
mins = header['axis mins']

labels = labels.astype(bool)
xx,yy,zz = labels.shape
slices = []

for z_slice in range(labels.shape[2]):
    current_slice = labels[:, :, z_slice]
    zoomed = zoom(current_slice, (0.5, 0.5))
    print("slice shape: " + str(zoomed.shape))
    slices.append(zoomed)

code.interact(local=locals())
labels = np.asarray(slices)
print("New shape: " + str(labels.shape))

skels = kimimaro.skeletonize(labels,
    teasar_params={
        'scale': 4,
        'const': 500, # physical units
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

"""
file = open("dendrite01_part2", "wb")
p.dump(allSkels, file)
"""
code.interact(local=locals())
for bb,ss in enumerate(skels):
    for i, skel in enumerate(ss.components()):
        x = skel.transform[0][3]
        y = skel.transform[1][3]
        filename = str(skel.id) + '-' + str(bb) +"_" + str(i)
        with open(filename + '.swc', 'wt') as f:
            f.write(skel.to_swc())
        swc = pySWC.Swc(filename + '.swc')
        swc.translate(x, y, 0)
        swc.scale(0.011)
        swc.save_file('c06_pt2_' + filename + '_trans.swc')