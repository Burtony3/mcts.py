# FOR REF
# Epoch = list of length m
# NAIF  = list of length n
import numpy as np
import const

def simulate(id, finalBody, maxNumFlyby, dvBudget, fblist):
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

    # ______________________________________________________________________________________________________
    # Stochastic Searching for End Condition

    # define t0 for the G array (note: t0 can changed b/c the first layer under the root are is a TOF sel.)
    # I have a concern of how we are storing states -burto3
    t0 = lineage[-2].state

    # G0 Array Constants     (G0 needs to be a property of mcts object)
    gLength = g0.shape[1]
    gHeight = g0.shape[0]

    # Finding Last Planet/Epoch Pair
    planet0 = node[lineage[1]].state # node[id].sometingwong
    epoch0 = node[lineage[0]].state  # node[id].lineages epoch thing

    # Pulling Current Constraints
    numFlyby = (len(lineage)/2)-1
    dvAcc = node[id].dvAcc

    # Exploring randomly until constraints violoated
    while const.check(planet0, dvAcc, numFlyby):

        # Randomly Choosing Planet & Epoch Pair
        randEpoch  = np.random.randint(gLength+1, size=1)[0]
        randPlanet = np.random.randint(gHeight+1, size=1)[0]

        # Translating to notation used in PyKep Lambert
        planet1 = fblist[randPlanet]
        epoch1 =  t0 + G0   # <-- G0 should be a global var (not literally adding, but giving G0 an epoch value)

        # Performing Lambert Calculation
        lambert(planet0, epoch0, planet1, epoch1)  # check syntax of inputs
        # lambert returns vinfdep and vinfarr

        # Setup for next loop through
        planet0 = planet1
        epoch0 = epoch1

        # Storing new checks
        dvAcc += dv  # summation of the dataout list
        numflyby += 1

    if currentBody is finalBody:
        X = # TODO: Calculate budget usage percentage (eq. 6 | pg. 771 | @izzo)
    else:
        X = 0
    return X


