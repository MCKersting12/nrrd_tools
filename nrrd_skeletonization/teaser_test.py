import pickle
import pySWC

file = open('dendrite01_part2', 'rb')

allskels = pickle.load(file)

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
