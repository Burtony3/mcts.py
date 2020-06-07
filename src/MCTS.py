import node as nodeNew
import spiceypy as spk
import numpy as np

# TODO: Actually import bsp/pck/tls files
spk.furnsh(stuff)

class MCTS(object):
    def __init__(self, finalBody, launchWindow, maxIters = 10000, dvBudget = 10, maxNumFlyby = float("inf"), detail = 30, flybyBodies = ["2", "3", "4", "5"]):
        """
        MCTS-Based Trajectory Combinational Solver
        
        Required Inputs:
            finalBody    :: String of the NAIF ID of the Arrival Body
            launchWindow :: List of strings with earliest & latest launch date (In form "MMM DD, YYYY")
        Optional Inputs:
            maxIters    :: Number of Nodes the Algorithm Creates (Default: 10,000)
            dvBudget    :: Maximum allowable Δv for entire sequence in km/s (Default: 10 km/s)
            maxNumFlyby :: Maximum number of flybys allowed before reaching final body (Default: inf)
            detail      :: Number of indicies for launch window & flyby dates (Default: 30?)
            flybyBodies :: List of NAIF ID Strings of planets available to be flown by (Default: ["2", "3", "4", "5"])
            ? → maxTof      :: Largest allowable mission elapsed time (Default: )
        Ouputs:
            ? → ID   :: Leaf Node of Highest Value
            ? → node :: List of all node objects
        """

        # Creating Root Node
        node = []                              # Allocates List
        node.append(nodeNew('3', None, 0, 0))  # Root of all Nodes (nodeNew will replace node.py)

        # FIXME:
        bodyDict = dict()

        # Create TOF array (G)
        # TODO: Adjust periods for input ephemeris
        T = np.array([88.0, 224.7, 365.2, 687.0, 4331, 10747, 30589, 59800]) * 84600
        G0 = []
        for i in range(len(T)):
            for j in range(detail):
                G0.append( j*detail*T(i)/360 )
        G0 = np.reshape(G0, (len(T), detail))


        #epochs = [self.sequence[1] + x * (self.l * int(T[int(self.sequence[-1])])) / 360 for x in range(9)]

        # Running to Completion
        for _ in range(maxIters):
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

if __name__ == "__main__":
    run inputs
