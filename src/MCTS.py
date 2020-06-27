class MCTS(object):
    def __init__(
        self,
        finalBody,  # String of NAIF ID
        launchWindow,  # Beginning and End Date String
        maxIters=10000,  # Maximum number of calculation iterations
        dvBudget=10,  # Maximum correction fuel budget
        maxNumFlyby=float("inf"),  # Maximum number of flyby's in a sequence
        detail=30,  # Number of indicies in launch & epoch nodes
        flybyBodies=["2", "3", "4", "5"],  # List of flyby body possibilities
        debug=False,
        oldTreeStyle=True,
    ):  # Debug printing

        # ------------------------- Package Installation ------------------------- #
        self.spk      = __import__("spiceypy")
        self.np       = __import__("numpy")
        self.math     = __import__("math")
        self.tabulate = __import__("tabulate")
        self.time     = __import__("time")
        self.rnd      = __import__("random")
        self.pk       = __import__("pykep")
        # ------------------------------------------------------------------------ #

        # Calculation Start Time
        self.startTime = self.time.time()

        self.oldTreeStyle = oldTreeStyle

        self.detail = detail

        # Furnshing Kernels
        self.spk.furnsh("../data/spk/naif0009.tls")
        self.spk.furnsh("../data/spk/de438.bsp")

        # ============================= DUMMY VALUES ============================= #
        self.P = flybyBodies  # Planets Venus through Saturn
        self.t0 = self.np.linspace(
            self.spk.utc2et(launchWindow[0]),
            self.spk.utc2et(launchWindow[1]),
            detail
        )
        # ======================================================================== # 

        # Initializing Tree
        self.node = []  # Allocating List
        if oldTreeStyle:
            self.node.append(nodeObj("3", parent=None, layer=1))  # Tree Root (id = 0)
            self.node[0].children = []
            for i in range(len(self.t0)):  # Looping through launch dates
                self.node.append(nodeObj(self.t0[i], layer=2))
                self.node[0].children.append(len(self.node) - 1)
        else:
            self.node.append(nodeObj(None, parent=None, layer=1))
            self.node[0].children = []
            for i in range(len(self.t0)):
                self.node.append(nodeObj(("3", self.t0[i]), layer=2))
                self.node[0].children.append(len(self.node) - 1)

        # Initializing Constraint Class
        const = constObj(finalBody, dvBudget, maxNumFlyby)

        # Converts NAIF ID's to row indexes for G0 array
        self.bodyDict = {}
        for i in range(len(self.P)):
            self.bodyDict[self.P[i]] = i

        # Create TOF array (G)
        # T = np.array([88.0, 224.7, 365.2, 687.0, 4331, 10747, 30589, 59800]) * 84600
        T = self.np.array([224.7, 365.2, 687.0, 4331]) * 84600
        Δ = 360/detail
        G0 = []
        for i in range(len(T)):
            for j in range(detail):
                G0.append((j+1)*Δ*T[i]/360)
        self.G0 = self.np.reshape(G0, (len(T), detail))


        # epochs = [self.sequence[1] + x * (self.l * int(T[int(self.sequence[-1])])) / 360 for x in range(9)]

        # Running to Completion
        for _ in range(maxIters):
            # Select
            id = self.select()  # Loops until finds leaf node & boolean if expand
            # print("Chosen node id = ", id, " | number of node visits = ")

            # Expand?
            if not self.node[id].children and (self.node[id].n != 0 or self.node[id].layer == 2):  # Expand if no children
                id = self.expand(id)

            
            self.simulate(id)  # <---- TEMP

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

    ## ———————————————————————————————————————————————————————————————————————————— ##
    ## —————————————————————————————— SUB-ROUTINES ———————————————————————————————— ##
    ## ———————————————————————————————————————————————————————————————————————————— ##

    def select(self):
        # Constant Initialization
        cp = 1 / self.math.sqrt(2.0)  # Cost Adjustment Parameter
        id = 0  # Starting at top of tree
        node = self.node  # Reassigning node for simplicity

        while node[id].children != None:
            # Updating Visits
            self.node[id].n += 1

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
                    print(ucb1, id, X, N, n, len(self.node))
                    ucb1.append(X + cp * self.math.sqrt(self.math.log1p(N) / n))

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
        self.node[id].n += 1

        np = self.np

        # Getting lineage to current node
        lineage = self.getLineage(id)

        # Shortenning node variable
        node = self.node

        # Changing expand node inputs based on node type
        if self.oldTreeStyle:
            if len(lineage) % 2 == 0:
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                # ~~~~~~~~~~~~~~ Planet Node ~~~~~~~~~~~~~~ #
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                states = self.P  # ← PLACEHOLDER
                dvAcc_ = None
                vInfIn_ = None
                vInfOut_ = None

            else:
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                # ~~~~~~~~~~~~~~~ Epoch Node ~~~~~~~~~~~~~~ #
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                lineage = self.getLineage(id)

                """
                # FIXME: Needs a dictionary in MCTS.py (bodyDict) to conv. bdy. str->num
                # Get row index of the G array from MCTS.py
                G0row = bodyDict(node[parent].state)
                # Get column that has the G state values (TOF vals.)
                states = G0[G0row]

                # TODO: Pruning nodes based on Periods and TOF
                # TOF eqn in page 771 Section 3
                # Should they be periods or ET?
                Ti_lower = 0.1*(node[lineage[3]].state + node[lineage[1]].state)
                Ti_upper = 2.0*(grandparent.state + parent.state)

                outputStates[:] = [states[i] if states[i] < Ti_upper and states[i] > Ti_lower for i in len(states)]
                """
                idx = self.bodyDict[self.node[id].state]
                states = self.t0[idx] + self.G0[idx, :]

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
        else:
            lineage = self.getLineage(id)
            states = []
            for i in range(len(self.P)):
                idx = self.bodyDict[self.P[i]]
                for j in range(len(self.t0)):
                    states.append((
                        self.P[i], 
                        self.node[lineage[-2]].state[1] + self.G0[idx, j],
                    ))

        # Getting Current Length of Node
        len_ = len(node)

        for i in range(len(states)):
            self.node.append(nodeObj(states[i], parent=lineage[0], layer=len(lineage) + 1))
            

            if i == 0:
                id = len(node)

        # Assigning new children ID's to parent node
        self.node[lineage[0]].children = [i for i in range(len_, len(node))]

        return id

    def simulate(self, id):
        lineage = self.getLineage(id)
        parent = lineage[1]
        self.node[id].n += 1
        self.node[id].cost = self.rnd.uniform(0, 1)
        
        # Departure Body State
        depState = self.spk.spkezr(
            self.node[parent].state[0],
            self.node[parent].state[1],
            'J2000',
            'NONE',
            '0'
        )[0]

        # Arrival Body State
        arrState = self.spk.spkezr(
            self.node[id].state[0],
            self.node[id].state[1],
            'J2000',
            'NONE',
            '0'
        )[0]

        # Calculating Time of Flight
        tof = self.node[id].state[1] - self.node[parent].state[1]


        if tof > 0:

            # Running Lambert Problem
            # print(depState[0:2],arrState[0:2],tof)
            l = self.pk.lambert_problem(
                r1 = depState[0:3],
                r2 = arrState[0:3],
                tof = tof,
                mu = self.pk.MU_SUN*1e-9
            )

            # Retrieving Lambert Calculation Velocities
            v0 = l.get_v1()
            v = l.get_v2()

            # Calculating V∞'s
            self.node[id].vInfOut = v0 - depState[3:6]
            self.node[id].vInfIn = v - arrState[3:6]

        

        

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
        time_ = self.time.time() - self.startTime
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
                children = tree[i].children
                if tree[i].children and len(tree[i].children) > 4:
                    children[4:] = ["..."]
                if tree[i].cost:
                    tree[i].cost = round(tree[i].cost, ndigits=3)
                cells.append(
                    [
                        i,
                        tree[i].layer,
                        tree[i].state,
                        # type_,
                        tree[i].n,
                        tree[i].cost,
                        tree[i].parent,
                        children,
                    ]
                )
            print(
                self.tabulate.tabulate(
                    cells,
                    ["id", "L", "State", "n", "Cost", "Par", "Children"]
                )
            )
        print(" ")
        print("EPOCH DEFINITION ARRAY (G0):\n", self.np.round(self.G0/84600, 3), "\n")
        print("RUN TIME:", round(time_, ndigits=5), " SECONDS")
        nLayer = []
        for i in range(len(tree)):
            if tree[i].layer > len(nLayer):
                nLayer.append(tree[i].n)
            else:
                nLayer[tree[i].layer - 1] += tree[i].n
        print("LAYER VISITS:", nLayer)
        print("EQUIVALENT RESOLUTION PARAMETER:", 360/self.detail, "\bᵒ")
        print(" ")


### ========================================================================== ###
### =============================== NODE CLASS =============================== ###
### ========================================================================== ###


class nodeObj:
    def __init__(self, state, parent=0, dvAcc=0, vInfIn=float("NaN"), layer=None):
        self.n = 0  # Visits in UCB1 Function
        # Estimated Reward in UCB1 Function (if -1 = dead node?)
        self.cost = []
        self.dvAcc = dvAcc  # Δv used up to that point
        # Establishes Parent above (Single Scalar of IDs)
        self.parent = parent
        self.children = None  # Children Below (Array of IDs)
        self.state = state  # Either NAIF ID String or Epoch
        self.vInfIn = vInfIn  # V∞ into CURRENT planet/epoch pair
        self.vInfOut = None  # V∞ out of PREVIOUS planet/epoch pair
        self.layer = layer

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
    mcts = MCTS(
        "5", ["Apr 01, 2020", "Jun 01, 2020"], maxIters=100, detail=8, debug=True, oldTreeStyle=False
    )

