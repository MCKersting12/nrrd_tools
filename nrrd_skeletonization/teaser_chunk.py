import nrrd
import kimimaro
import pickle as p

labels, header = nrrd.read('VCN_c06_dendrite01_part1.nrrd')
print(labels.shape)
mins = header['axis mins']

labels = labels.astype(bool)

xx,yy,zz = labels.shape
halfx = xx/2
halfy = yy/2

allSkels = []
for x in range(2):
  for y in range(2):
    block = labels[int(x*halfx):int(x*halfx + halfx), int(y*halfy):int(y*halfy + halfy), :]
    print(block.shape)
    skels = kimimaro.skeletonize(block,
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
    if not skels:
      continue
    for skel in skels.values():
      skel.transform[0][3] = int(x*halfx) + mins[0]
      skel.transform[1][3] = int(y*halfy) + mins[1]
      skel.transform[2][3] = mins[2]

    allSkels.append(skels)

file = open("dendrite01_part2", "wb")
p.dump(allSkels, file)
"""
newSkels = []

for ii,skeldict in enumerate(allskels):
    for jj,comp in enumerate(skeldict):
        a = skeldict[comp]
        x = a.transform[0][3]
        y = a.transform[1][3]
        newSkels.append(a)

for bb,ss in enumerate(newSkels):
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
"""
