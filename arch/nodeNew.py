# Class Creation
class node:

    # Blank Node Creation
    def __init__(self, state, parent = 0, dvAcc = 0, vInfIn = float("NaN"), vInfOut = float("NaN")):
        self.n = 0              # Visits in UCB1 Function
        self.cost = None        # Estimated Reward in UCB1 Function (if -1 = dead node?)
        self.dvAcc = dvAcc      # Δv used up to that point
        self.parent = parent    # Establishes Parent above (Single Scalar of IDs)
        self.child = None       # Children Below (Array of IDs)
        self.state = state      # Either NAIF ID String or Epoch
        self.vInfIn = vInfIn    # V∞ into CURRENT planet/epoch pair
        self.vInfOut = vInfOut  # V∞ out of PREVIOUS planet/epoch pair
