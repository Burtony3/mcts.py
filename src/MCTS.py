class MCTS(object):
    def __init__(
        self,
        finalBody,  # String of NAIF ID
        launchWindow,  # Beginning and End Date String
        maxIters=10000,  # Maximum number of calculation iterations
        ΔvBudget=10,  # Maximum correction fuel budget
        maxNumFlyby=float("inf"),  # Maximum number of flyby's in a sequence
        detail=30,  # Number of indicies in launch & epoch nodes
        flybyBodies=['2', '3', '4'],  # List of flyby body possibilities
        debug=False,
        oldTreeStyle=True,
        frame = 'ECLIPJ2000',
        abcorr = 'NONE',
        cntrBody = '0',
        C3max = 100
    ):  # Debug printing

        # ------------------------- PACKAGE INSTALLATION ------------------------- #
        self.spk      = __import__("spiceypy")
        self.np       = __import__("numpy")
        self.math     = __import__("math")
        self.tabulate = __import__("tabulate")
        self.time     = __import__("time")
        self.rnd      = __import__("random")
        self.pk       = __import__("pykep")
        # ------------------------------------------------------------------------ #

        # TAKING START TIME
        self.startTime = self.time.time()

        # SAVING INPUTS TO OBJECT
        self.finalBody = finalBody
        self.oldTreeStyle = oldTreeStyle
        self.detail = detail
        self.frame = frame
        self.abcorr = abcorr
        self.cntrBody = cntrBody
        self.C3max = C3max

        # FURNISHING KERNELS
        self.spk.furnsh("../data/spk/naif0009.tls")
        self.spk.furnsh("../data/spk/de438.bsp")

        # PARSING USER INPUTS

        self.P = flybyBodies + [finalBody]  # Planets Venus through Saturn
        # ============================= DUMMY VALUES ============================= #
        self.t0 = self.np.linspace(
            self.spk.utc2et(launchWindow[0]),
            self.spk.utc2et(launchWindow[1]),
            detail
        )
        # ======================================================================== # 

        # INITIALIZING DICTIONARIES
        self.pkP = {'1': 'mercury', '2': 'venus',   '3': 'earth', '4': 'mars', '5': 'jupiter', '6': 'saturn',  '7': 'uranus',  '8': 'neptune'}
        self.tau = {'1': 87.97,     '2': 224.7,     '3': 365.25,  '4': 687,    '5': 4331,      '6': 10747,     '7': 30589,     '8': 59800}
        self.col = {'1': '#5a6165', '2': '#955c1c', '3': 'b',     '4': 'r',    '5': '#9e7a5b', '6': '#c9c0a9', '7': '#8eb0b8', '8': '#4d80d7'}

        # INITIALIZING TREE
        self.treeArt()
        self.node = []  # Allocating List
        if oldTreeStyle:
            self.node.append(nodeObj("3", parent=None, layer=1))  # Tree Root (id = 0)
            self.node[0].children = []
            for i in range(len(self.t0)):  # Looping through launch dates
                self.node.append(nodeObj(self.t0[i], layer=2, Δv = ΔvBudget))
                self.node[0].children.append(len(self.node) - 1)
        else:
            self.node.append(nodeObj(None, parent=None, layer=1))     # Creating Root Node
            self.node[0].children = []
            for i in range(len(self.t0)):
                self.node.append(nodeObj(('3', self.t0[i]), layer=2, Δv = ΔvBudget))
                self.node[0].children.append(len(self.node) - 1)

        # Initializing Constraint Class
        # const = constObj(finalBody, dvBudget, maxNumFlyby)

        # Converts NAIF ID's to row indexes for G0 array
        self.bodyDict = {}
        for i in range(len(self.P)):
            self.bodyDict[self.P[i]] = i

        # Create TOF array (G)
        # T = np.array([88.0, 224.7, 365.2, 687.0, 4331, 10747, 30589, 59800]) * 84600
        T = [self.tau[self.P[i]] for i in range(len(self.P))]
        Δ = 360/detail
        G0 = []
        for i in range(len(T)):
            for j in range(detail):
                G0.append(86400*(j+1)*Δ*T[i]/360)
        self.G0 = self.np.reshape(G0, (len(T), detail))


        # epochs = [self.sequence[1] + x * (self.l * int(T[int(self.sequence[-1])])) / 360 for x in range(9)]

        # Running to Completion
        for _ in range(maxIters):
            # Select
            id = self.select()  # Loops until finds leaf node & boolean if expand
            # print("Chosen node id = ", id, " | number of node visits = ")

            # Expand?
            if not self.node[id].children and (self.node[id].n != 0): #or self.node[id].layer == 2):  # Expand if no children
                id = self.expand(id)

            
            self.simulate(id)  # <---- TEMP

        """
            # Simulate
            if not string(node[id].state):
                value = simulate(id, finalBody, maxNumFlyby, dvBudget, fblist)     #update function call's inputs

                # Backprop
                backprop(id, value)
        """

        self.spk.kclear()

        ### DEBUGG PRINTING ###
        if debug:
            self.debug()

    ## ———————————————————————————————————————————————————————————————————————————— ##
    ## —————————————————————————————— SUB-ROUTINES ———————————————————————————————— ##
    ## ———————————————————————————————————————————————————————————————————————————— ##

    def select(self):
        # CONSTANT INITIALIZATION
        cp = 1 / self.math.sqrt(2.0)  # Cost Adjustment Parameter
        id = 0                        # Starting at top of tree
        node = self.node              # Reassigning node for clarity

        # PATHING TO HIGHEST VALUE LEAF
        while node[id].children != None:
            # UPDATING VISITS
            self.node[id].n += 1

            # INITIALIZING USEFUL VARIABLES
            ucb1 = []
            len_ = len(node[id].children)


            # LOOPING THROUGH AVAILABLE CHILDREN
            for i in range(len_):
                X = node[node[id].children[i]].X
                N = node[id].n
                n = node[node[id].children[i]].n
                isTerminal = node[node[id].children[i]].isTerminal

                if n == 0 and not isTerminal:  # Selecting first unvisited node
                    id = node[id].children[i]
                    break
                elif not isTerminal:           # Calculating Value
                    # print(ucb1, id, X, N, n, len(self.node))
                    # self.debug()
                    ucb1.append(X + cp * self.math.sqrt(self.math.log1p(N) / n))
                else:                          # Skipping Terminal Node
                    ucb1.append(0)

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

        # Getting lineage to current node
        lineage = self.getLineage(id)

        # Shortenning node variable
        node = self.node

        lineage = self.getLineage(id)
        states = []
        for i in range(len(self.P)):
            idx = self.bodyDict[self.P[i]]
            # ← Node Pruning
            for j in range(len(self.t0)):
                states.append((
                    self.P[i], 
                    self.node[lineage[-2]].state[1] + self.G0[idx, j],
                ))

        # Getting Current Length of Node
        len_ = len(node)

        for i in range(len(states)):
            self.node.append(nodeObj(
                states[i], 
                parent=lineage[0], 
                layer=len(lineage) + 1, 
                isTerminal= states[i][0] == self.finalBody,
                Δv = self.node[lineage[1]].Δv))

            if self.node[-1].layer == 3:
                self.debug()
                l0 = self.computeLambert(len_+ i)
                Δv = self.computeΔv(len_+ i, l0)
                self.node[-1].Δv -= Δv

            if i == 0:
                id = len(node)

        # Assigning new children ID's to parent node
        self.node[lineage[0]].children = [i for i in range(len_, len(node))]

        return id

    def simulate(self, id):
        lineage = self.getLineage(id)
        parent = lineage[1]
        self.node[id].n += 1
        self.node[id].X = self.rnd.uniform(0, 1)
        
        # Departure Body State
        """
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
        """
        

        

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

    def computeLambert(self, id):
        # IMPORTING TOOLBOXES
        spk = self.spk
        pk = self.pk

        # GETTING PARENT
        parent = self.node[id].parent

        # GETTING STATES
        p0, et0 = self.node[parent].state
        p1, et1 = self.node[id].state

        # GETTING SPICE STATE
        s0 = spk.spkezr(p0, et0, self.frame, self.abcorr, self.cntrBody)[0]
        s1 = spk.spkezr(p1, et1, self.frame, self.abcorr, self.cntrBody)[0]

        # COMPUTING LAMBERT
        l = pk.lambert_problem(
            r1 = s0[0:3]*1000,
            r2 = s1[0:3]*1000,
            tof = et1 - et0,
            mu = pk.MU_SUN
        )

        return l

    def computeΔv(self, id, l0, l1 = None):
        # l₀ :: Lambert Arc from Parent → ID
        # l₁ ::                  Grandparent → Parent
        parent = self.node[id].parent
        vp = self.spk.spkezr(
            self.node[parent].state[0],
            self.node[parent].state[1],
            self.frame,
            self.abcorr,
            self.cntrBody
        )[0][3:]

        # CALCULATING Δv
        if not l1: # Returns Launch C3
            v0 = self.np.array(l0.get_v1()[0])/1000
            vinf = v0 - vp
            Δv = max(0, self.np.linalg.norm(vinf) - self.math.sqrt(self.C3max))
        else:
            vi = self.np.array(l1.get_v2()[0])/1000 - vp
            vo = self.np.array(l0.get_v1()[0])/1000 - vp
            Δv = self.np.linalg.norm(vo - vi)
        return Δv


    def treeArt(self):
        print("\n\n        . . .")
        print("      .        .  .     ..    .")
        print("   .                 .         .  .")
        print("                   .")
        print("                  .                ..")
        print("  .          .            .              .")
        print("  .            '.,        .               .   Algorithm Developed for use in AAS 20-686:")
        print("  .              'b      *                    Implementation of Deep Space Maneuvers in ")
        print("   .              '$    #.                ..  Broad Search Trajectories Using Monte Carlo")
        print("  .    .           $:   #:               .    Tree Search")
        print("..      .  ..      *#  @):        .   . .")
        print("             .     :@,@):   ,.**:'   .")
        print(" .      .,         :@@*: ..**'      .   .     Created by: Burton Yale, Rohan Patel,")
        print("          '#o.    .:(@'.@*\"\'  .                           & Jehosafat Cabrera")
        print("  .  .       'bq,..:,@@*'   ,*      .  .")
        print("             ,p$q8,:@)'  .p*'      .")
        print("      .     '  . '@@Pp@@*'    .  .")
        print("       .  . ..    Y7'.'     .  .")
        print("                 :@):.")
        print("                .:@:'.")
        print("______________.::(@:.______________________ \n\n")

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
                if tree[i].X:
                    tree[i].X = round(tree[i].X, ndigits=3)
                cells.append(
                    [
                        i,
                        tree[i].layer,
                        tree[i].isTerminal,
                        tree[i].state,
                        # type_,
                        tree[i].n,
                        tree[i].X,
                        tree[i].Δv,
                        tree[i].parent,
                        # children,
                    ]
                )
            print(
                self.tabulate.tabulate(
                    cells,
                    ["id", "L", "Leaf", "State", "n", "X", "Δv", "Par"]
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
    def __init__(self, state, parent=0, Δv=0, vInfIn=float("NaN"), layer=None, isTerminal=False):
        self.n = 0  # Visits in UCB1 Function
        # Estimated Reward in UCB1 Function (if -1 = dead node?)
        self.X = []
        self.Δv = Δv  # Δv used up to that point
        # Establishes Parent above (Single Scalar of IDs)
        self.parent = parent
        self.children = None  # Children Below (Array of IDs)
        self.state = state  # Either NAIF ID String or Epoch
        self.layer = layer
        self.isTerminal = isTerminal

    def update(self):
        print("dif")


### ================================================================================ ###
### ============================== CONSTRAINTS CLASS =============================== ###
### ================================================================================ ###

"""
class constObj:
    def __init__(self, finalBody, ΔvBudget, maxNumFlyby):
        self.dv = dvBudget
        self.flyby = maxNumFlyby
        self.finalBody = finalBody

    def check(self, finalBody, dv, numFlyby):
        constList = [finalBody != self.finalBody, dv < self.dv, numFlyby < self.flyby]
        return all(constList)
"""

### ============================================================================== ###
### =============================== INITIALIZATION =============================== ###
### ============================================================================== ###

if __name__ == "__main__":
    mcts = MCTS(
        '5', ["Apr 01, 2020", "Jun 01, 2020"], maxIters=100, detail=8, debug=True, oldTreeStyle=False
    )

