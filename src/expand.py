import numpy as np

def expand(id):
    lineage = getLineage(id)
    if len(lineage) % 2 == 0:
        states = [1, 2, 3, 4, 5, 6, 7, 8]
        dvAcc_ = None
        vInfIn_ = None
        vInfOut_ = None
    else:
        # We have to code the G function here
        states = G function # Creating the list of states
        for i in range(len(states)):
        # Pruning nodes based on Periods and TOF
        grandparent = lineage[1]
        dvAcc_ = node[grandparent].dvAcc
        vInfIn_ = np.zeros([3 1])
        vInfOut_ = np.zeros([3 1])

    for i in range(len(states)):
        node.append(states[i], parent=id, dvAcc=dvAcc_, vInfIn = vInfIn_, vInfOut = vInfOut_)
