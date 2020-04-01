# FOR REF
# Epoch = list of length m
# NAIF  = list of length n
import numpy as np

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

    # Global G0 Array Constants     (G0 needs to be global!)
    gLength = g0.shape[1]
    gHeight = g0.shape[0]

    planet0 = # node[id].sometingwong
    epoch0 = # node[id].lineages epoch thing

    dataout = []
    numflyby = (len(lineage)/2)-1

    # Exploring randomly in Simulate
    while checkConstraints():

        randEpoch  = np.random.randint(gLength+1, size=1)[0]
        randPlanet = np.random.randint(gHeight+1, size=1)[0]

        planet1 = fblist[randPlanet]
        epoch1 =  t0 + G0   # <-- G0 should be a global var (not literally adding, but giving G0 an epoch value)

        #lambert(planet0, epoch0, planet1, epoch1)  # check syntax of inputs
        # lambert returns vinfdep and vinfarr

        planet0 = planet1
        epoch0 = epoch1

        dataout.append()    #<--- append DV from lambert here (no checks for max dv
        dvAcc = # summation of the dataout list

        numflyby += 1


