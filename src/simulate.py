# FOR REF
# Epoch = list of length m
# NAIF  = list of length n

def simulate(id, finalBody, maxNumFlyby, dvBudget):
    # Getting Current Node Lineage
    lineage = getLineage(id)

    # Running Current Pair Lambert
    vInfO, vInfI = lambert(node[lineage[3]].state,
                           node[lineage[2]].state,
                           node[lineage[1]].state,
                           node[lineage[0]].state)

    # Saving to Node
    node[id].vInfOut = vInfO
    node[id].vInfIn = vInfI

    # Δv's at previous body
    dV = np.linalg.norm(np.subtract(vInfO, node[lineage[2]].vInfIn))
    node[id].dvAcc += dV

    # Stochastic Searching for End Condition
    while curBody not finalBody and dvAcc < dvBudget and numFlyby < maxNumFlyby:
