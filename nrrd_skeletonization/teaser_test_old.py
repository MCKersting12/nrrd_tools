import nrrd
import code
import kimimaro
import scipy
import pickle
import pySWC

file = open('dendrite01_part2', 'rb')

# dump information to that file
allskels = pickle.load(file)

newSkels = []
uuid = 0
for ii,skeldict in enumerate(allskels):
    for jj,comp in enumerate(skeldict):
        #code.interact(local=locals())
        #print(skeldict[comp].vertices)

        a = skeldict[comp]
        x = a.transform[0][3]
        y = a.transform[1][3]

       # code.interact(local=locals())

        filename = str(uuid) + "_" + str(a.id) + '-' + str(ii) + "_" + str(jj)
        with open(filename + '.swc', 'wt') as f:
            f.write(a.to_swc())
        swc = pySWC.Swc(filename + '.swc')
        swc.translate(x, y, 0)
        swc.scale(0.011)
        swc.save_file(filename + "_trans.swc")

        newSkels.append(a)
        uuid += 1

for bb,ss in enumerate(newSkels):
    for i, skel in enumerate(ss.components()):
        filename = str(skel.id) + '-' + str(bb) +"_" + str(i) + '.swc'
        with open(filename, 'wt') as f:
            f.write(skel.to_swc())

"""
from cloudvolume import PrecomputedSkeleton

skel = PrecomputedSkeleton.simple_merge(newSkels).consolidate()
skel = kimimaro.postprocess(
  skel,
  dust_threshold=1000, # physical units
  tick_threshold=3500 # physical units
)

skeleton = skel

with open('theWholeSkel_2.swc', 'wt') as f:
    f.write(skeleton.to_swc())
"""

# Split input skeletons into connected components and
# then join the two nearest vertices within `radius` distance
# of each other until there is only a single connected component
# or no pairs of points nearer than `radius` exist.
# Fuse all remaining components into a single skeleton.
#skel = kimimaro.join_close_components([skel1, skel2], radius=1500) # 1500 units threshold

#skel = kimimaro.join_close_components([skel1, skel2], radius=None) # no threshold

#for i, skel in enumerate(skeleton.components()):






