from node import * as nodeNew
import spiceypy as spk

spk.furnsh(stuff)

class MCTS(object):
    def __init__(maxIter, dvBudget, maxNumFlyby, launchWindow, finalBody; kwargs...):
        # Input Parse


        # Creating Root Node
        node = []                              # Allocates List
        node.append(nodeNew('3', None, 0, 0))  # Root of all Nodes (nodeNew will replace node.py)

        # Running to Completion
        for i in range(maxIters):
            # Select
            id, expandBool = select()   # Loops until finds leaf node & boolean if expand

            # Expand?
            if expandBool is True:
                id = expand(id, finalBody)

            # Simulate
            if not string(node[id].state):
                value = simulate(id, finalBody, maxNumFlyby, dvBudget)

                # Backprop
                backprop(id, value)

        # Returning Results
        return node

spk.kclear()
