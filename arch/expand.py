import numpy as np

def expand(id):
    lineage = getLineage(id)
    if len(lineage) % 2 == 0:
        states = flybyBodies
        states.append(finalBody)
        dvAcc_ = None
        vInfIn_ = None
        vInfOut_ = None
    else:

        # This condition is triggered if we are expanding from a planet node. Need to expand to a TOF node

        parent = lineage[0]         # Parent is a Planet Node (that expand is creating children for)
        grandparent = lineage[2]    # Also a Planet Node (used for Tau_i-1)


        # FIXME: Needs a dictionary in MCTS.py (bodyDict) to conv. bdy. str->num
        G0row = bodyDict(node[parent].state)           # Get row index of the G array from MCTS.py
        states = G0[G0row]                             # Get column that has the G state values (TOF vals.)


        # TODO: Pruning nodes based on Periods and TOF
        # TOF eqn in page 771 Section 3
        Ti_lower = 0.1*(grandparent.state + parent.state)
        Ti_upper = 2.0*(grandparent.state + parent.state)

        #outputStates[:] = [states[i] if states[i] < Ti_upper and states[i] > Ti_lower for i in len(states)]

        outputStates = []
        for i in range(len(states)):
            if states[i] > Ti_lower && states[i] < Ti_upper
                outputStates.append(states[i])

        outputStates.append(Ti_upper)
        outputStates.insert(0, Ti_lower)


        dvAcc_ = node[lineage[1]].dvAcc
        vInfIn_ = np.zeros([3 1])
        vInfOut_ = np.zeros([3 1])


    for i in range(len(states)):

        node.append(outputStates[i], parent=id, dvAcc=dvAcc_, vInfIn = vInfIn_, vInfOut = vInfOut_)

        if i == 0:
            id = len(node)

    return id