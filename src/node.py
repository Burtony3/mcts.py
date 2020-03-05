# Importing
import numpy as np
from sequenceEnv import * as env

# Class Creation
class node(object):

    # Blank Node Creation
    def __init__(self, state, parent=0, dvAcc = 0, vInfIn = np.zeros([3 1]), vInfOut):
        self.n = 0              # Visits in UCB1 Function
        self.cost = None        # Estimated Reward in UCB1 Function (if -1 = dead node)
        self.dvAcc = dvAcc      # Δv used up to that point
        self.parent = parent    # Establishes Parent above (Single Scalar of IDs)
        self.child = None       # Children Below (Array of IDs)
        self.state = state      # Either NAIF ID String or Epoch
        self.isLeaf = True      # Is end of branch
        self.vInfIn = vInfIn    # V∞ into CURRET planet/epoch pair
        self.vInfOut = vInfOut  # V∞ out of PREVIOUS planet/epoch pair
