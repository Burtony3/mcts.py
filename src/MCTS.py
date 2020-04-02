from node import * as nodeNew
import spiceypy as spk
import numpy as np

# TODO: Actually import bsp/pck/tls files
spk.furnsh(stuff)

class MCTS(object):
    def __init__(maxIter, dvBudget, maxNumFlyby, launchWindow, finalBody, glength; kwargs...):
        # Input Parse
        # TODO: Make optional inputs

        # Flyby Body List
        fblist = ['1','2','3','4','5','6','7','8','9']

        # Creating Root Node
        node = []                              # Allocates List
        node.append(nodeNew('3', None, 0, 0))  # Root of all Nodes (nodeNew will replace node.py)

        # Create TOF array (G)
        # TODO: Adjust periods for input ephemeris
        T = np.array([88.0, 224.7, 365.2, 687.0, 4331, 10747, 30589, 59800]) * 84600
        G0 = []
        # the height of the G0 array must be the same value as the fblist (by def.)


        #epochs = [self.sequence[1] + x * (self.l * int(T[int(self.sequence[-1])])) / 360 for x in range(9)]

        # Running to Completion
        for i in range(maxIters):
            # Select
            id, expandBool = select()   # Loops until finds leaf node & boolean if expand

            # Expand?
            if expandBool is True:
                id = expand(id, finalBody)

            # Simulate
            if not string(node[id].state):
                value = simulate(id, finalBody, maxNumFlyby, dvBudget, fblist)     #update function call's inputs

                # Backprop
                backprop(id, value)

        # Returning Results
        return node

spk.kclear()
