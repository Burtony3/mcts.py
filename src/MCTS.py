class MCTS(object):
    def __init__(self, 
        finalBody,                          # String of NAIF ID
        launchWindow,                       # Beginning and End Date String
        maxIters    = 10000,                # Maximum number of calculation iterations
        dvBudget    = 10,                   # Maximum correction fuel budget
        maxNumFlyby = float("inf"),         # Maximum number of flyby's in a sequence
        detail      = 30,                   # Number of indicies in launch & epoch nodes
        flybyBodies = ["2", "3", "4", "5"], # List of flyby body possibilities
        debug       = False):               # Debug printing

        # ------------------------- Package Installation ------------------------- # 
        self.spk = __import__('spiceypy')
        self.np = __import__('numpy')
        self.math = __import__('math')
        self.tabulate = __import__('tabulate')
        self.time = __import__('time')
        # ------------------------------------------------------------------------ #

        # Calculation Start Time
        self.startTime = self.time.time()

        # ========================= DUMMY VALUES ========================= #
        self.planets = ['2', '3', '4', '5', '6']            # Planets Venus through Saturn
        self.dates = [2459017.5 + i*30 for i in range(50)]  # jdates from 2020-06-16 to 2022-06-16
        # ================================================================ #

        # Initializing Tree
        self.node = []                                # Allocating List
        self.node.append(nodeObj('3', parent=None, layer=1))   # Tree Root (id = 0)
        self.node[0].children = []
        for i in range(len(self.dates)):              # Looping through launch dates
            self.node.append(nodeObj(self.dates[i], layer=2))
            self.node[0].children.append(len(self.node)-1)
            # len_ = len(self.node)-1                   # Getting id of launch date
            # for j in range(len(self.planets)):        # Looping through flyby planets
            #     self.node.append(nodeObj(self.planets[j], parent=len_, layer=3))
            # self.node[len_].children = [i for i in range(len_+1, len(self.node))]
            # self.node[0].children.append(len_)

        # Initializing Constraint Class
        const = constObj(finalBody, dvBudget, maxNumFlyby)

        """
        # TODO: Convert NAIF ID's to row indexes for G0 array
        bodyDict = dict()

        # Create TOF array (G)
        # TODO: Adjust periods for input ephemeris
        T = np.array([88.0, 224.7, 365.2, 687.0, 4331, 10747, 30589, 59800]) * 84600
        G0 = []
        for i in range(len(T)):
            for j in range(detail):
                G0.append( j*detail*T(i)/360 )
        G0 = np.reshape(G0, (len(T), detail))


        # epochs = [self.sequence[1] + x * (self.l * int(T[int(self.sequence[-1])])) / 360 for x in range(9)]
        """

        # Running to Completion
        for _ in range(maxIters):
            # Select
            id = self.select()   # Loops until finds leaf node & boolean if expand
            # print("Chosen node id = ", id, " | number of node visits = ")
            
            # Expand?
            if not self.node[id].children and self.node[id].n != 0: # Expand if no children
                id = self.expand(id)

            self.simulate(id)

        """
            # Simulate
            if not string(node[id].state):
                value = simulate(id, finalBody, maxNumFlyby, dvBudget, fblist)     #update function call's inputs

                # Backprop
                backprop(id, value)
        """

        ### DEBUGG PRINTING ###
        if debug:
            self.debug()

    ## ——————————————————————————————————————————————————————————————————————————— ##
    ## —————————————————————————————— SUB-ROUTINES ——————————————————————————————— ##
    ## ——————————————————————————————————————————————————————————————————————————— ##
        
    def select(self):
        # Constant Initialization
        cp = 1/self.math.sqrt(2.0)  # Cost Adjustment Parameter
        id = 0                 # Starting at top of tree
        node = self.node       # Reassigning node for simplicity

        while node[id].children != None:
            # Updating Visits
            self.node[id].n = self.node[id].n + 1

            ucb1 = []
            len_ = len(node[id].children)
            # Checks Node's children
            for i in range(len_):
                X = node[node[id].children[i]].cost
                N = node[id].n
                n = node[node[id].children[i]].n

                if n == 0:  # Immediatly selects first unvisited node
                    id = node[id].children[i]
                    break
                else:
                    ucb1.append(X + cp*self.math.sqrt(self.math.log1p(N)/n))

            # Checking whether UCB1 is empty or not
            if len(ucb1) is len_:
                # Find UCB w/ max value
                indexChild = ucb1.index(max(ucb1))

                # Set the id to this child's index in the Node list
                id = node[id].children[indexChild]


        """
        # Skips expand step if already has visits
        # TODO: Check if needed?
        if node[id].n == 0:
            expandBool = False
        else:
            expandBool = True
        """

        return id

    def expand(self, id):
        np = self.np

        # Getting lineage to current node
        lineage = self.getLineage(id)

        # Shortenning node variable
        node = self.node
        
        # Changing expand node inputs based on node type
        if len(lineage) % 2 == 0:
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # ~~~~~~~~~~~~~~ Planet Node ~~~~~~~~~~~~~~ #
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

            states = self.planets # ← PLACEHOLDER
            dvAcc_ = None
            vInfIn_ = None
            vInfOut_ = None

        else:
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # ~~~~~~~~~~~~~~~ Epoch Node ~~~~~~~~~~~~~~ #
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

            """
            # FIXME: Needs a dictionary in MCTS.py (bodyDict) to conv. bdy. str->num
            # Get row index of the G array from MCTS.py
            G0row = bodyDict(node[parent].state)
            # Get column that has the G state values (TOF vals.)
            states = G0[G0row]
            """
            
            ### 🚩 STOPPED HERE 🚩

            """
            # TODO: Pruning nodes based on Periods and TOF
            # TOF eqn in page 771 Section 3
            # Should they be periods or ET?
            Ti_lower = 0.1*(node[lineage[3]].state + node[lineage[1]].state)
            Ti_upper = 2.0*(grandparent.state + parent.state)

            outputStates[:] = [states[i] if states[i] < Ti_upper and states[i] > Ti_lower for i in len(states)]
            """

            states = self.dates # ← Placeholder

            """
            # Checks for if outputStates meets Ti bounds, then adds it to states
            for i in range(len(states)):
                if (states[i] > Ti_lower) and (states[i] < Ti_upper):
                    outputStates.append(states[i])

            # Inserting bounds as possible dates
            outputStates.append(Ti_upper)
            outputStates.insert(0, Ti_lower)
            """

            dvAcc_ = node[lineage[1]].dvAcc
            vInfIn_ = np.zeros([3, 1])
            vInfOut_ = np.zeros([3, 1])

        # Getting Current Length of Node
        len_ = len(node)

        for i in range(len(states)):
            self.node.append(nodeObj(states[i], parent=lineage[0], layer=len(lineage)+1))

            if i == 0:
                id = len(node)

        # Assigning new children ID's to parent node
        self.node[lineage[0]].children = [i for i in range(len_, len(node))]

        return id

    def simulate(self, id):
        self.node[id].n = self.node[id].n + 1

    ## -------------------------------------------------------------------------------- ##
    ## ------------------------------- HELPER FUNCTIONS ------------------------------- ##
    ## -------------------------------------------------------------------------------- ##

    def getLineage(self, id):
        lineage = []
        while id != 0:
            lineage.append(id)
            id = self.node[id].parent
        lineage.append(0)

        # Returns most recent node @ idx = 0 & oldest parent @ idx = -1
        return lineage

    def debug(self):
        print(" ")
        print("### =================================================== ###")
        print("### ==================== DEBUGGING ==================== ###")
        print("### =================================================== ###")
        print(" ")
        cells = []
        tree = self.node
        if len(tree) > 500:
            r = input("TREE EXCEEDS 500 NODES CONTINUE? (Y/N): ")
            r = r == "Y" or r == "y"
        else:
            r = True
        if r:
            print("↓TREE PRODUCED↓")
            for i in range(len(tree)):
                if type(tree[i].state) is str:
                    type_ = "Planet"
                else:
                    type_ = "Epoch"
                cells.append([i, 
                    tree[i].layer, 
                    tree[i].state, 
                    type_, 
                    tree[i].n, 
                    tree[i].cost, 
                    tree[i].parent, 
                    tree[i].children])
            print(self.tabulate.tabulate(cells, ["id", "Layer", "State", "Type", "n", "Cost", "Par", "Children"]))
        print(" ")
        print("RUN TIME: ", round(self.time.time() - self.startTime, ndigits=5), " SECONDS")
        # print("MAX LAYER: ", max(layer for tree in tree))
        # print("ITERATIONS: ", maxIters)
        print(" ")



### ========================================================================== ###
### =============================== NODE CLASS =============================== ###
### ========================================================================== ###

class nodeObj:
    def __init__(self, state, parent=0, dvAcc=0, vInfIn=float("NaN"), layer = None):
        rnd = __import__('random')
        self.n = 0              # Visits in UCB1 Function
        # Estimated Reward in UCB1 Function (if -1 = dead node?)
        self.cost = rnd.uniform(0, 1)
        self.dvAcc = dvAcc      # Δv used up to that point
        # Establishes Parent above (Single Scalar of IDs)
        self.parent = parent
        self.children = None    # Children Below (Array of IDs)
        self.state = state      # Either NAIF ID String or Epoch
        self.vInfIn = vInfIn    # V∞ into CURRENT planet/epoch pair
        self.vInfOut = None  # V∞ out of PREVIOUS planet/epoch pair
        self.layer = layer

    def visit(self):
        self.n = self.n+1

    def update(self):
        print("dif")

### ================================================================================ ###
### ============================== CONSTRAINTS CLASS =============================== ###
### ================================================================================ ###

class constObj:
    def __init__(self, finalBody, dvBudget, maxNumFlyby):
        self.dv = dvBudget
        self.flyby = maxNumFlyby
        self.finalBody = finalBody

    def check(self, finalBody, dv, numFlyby):
        constList = [finalBody != self.finalBody, dv < self.dv, numFlyby < self.flyby]
        return all(constList)

### ============================================================================== ###
### =============================== INITIALIZATION =============================== ###
### ============================================================================== ###

if __name__ == "__main__":
    # TODO: Actually import bsp/pck/tls files
    # spk.furnsh(stuff)


    mcts = MCTS(
        '5', 
        ['Apr 01, 2020', 'Jun 01, 2020'], 
        maxIters=500, 
        debug = True)
