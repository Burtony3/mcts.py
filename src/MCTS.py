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
        self.plt      = __import__("matplotlib")
        self.plt      = self.plt.pyplot
        # ------------------------------------------------------------------------ #

        # TAKING START TIME
        self.startTime = self.time.time()

        # SAVING INPUTS TO OBJECT
        self.finalBody = finalBody
        self.flybyBodies = flybyBodies
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
            self.node.append(nodeObj("3", parent=None, layer=1, Δv = ΔvBudget))  # Tree Root (id = 0)
            self.node[0].children = []
            for i in range(len(self.t0)):  # Looping through launch dates
                self.node.append(nodeObj(self.t0[i], layer=2, Δv = ΔvBudget))
                self.node[0].children.append(len(self.node) - 1)
        else:
            self.node.append(nodeObj(None, parent=None, layer=1, Δv = ΔvBudget))     # Creating Root Node
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

        ### DEBUGG PRINTING ###
        if debug:
            self.debug()

        # self.spk.kclear()

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
                # print(id)
                X = node[node[id].children[i]].X
                N = node[id].n
                n = node[node[id].children[i]].n
                isTerminal = node[node[id].children[i]].isTerminal

                if n == 0 and not isTerminal:  # Selecting first unvisited node
                    id = node[id].children[i]
                    break
                elif isTerminal:           # Calculating Value
                    # print(ucb1, id, X, N, n, len(self.node))
                    # self.debug())
                    ucb1.append(0)
                else:                          # Skipping Terminal Node
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

        # Getting lineage to current node
        lineage = self.getLineage(id)

        # Shortenning node variable
        node = self.node

        lineage = self.getLineage(id)
        states = []
        P_ = self.flybyBodies + [self.finalBody]
        P_.remove(self.node[id].state[0])
        for i in range(len(P_)):
            idx = self.bodyDict[P_[i]]
            # ← Node Pruning
            for j in range(len(self.t0)):
                states.append((
                    P_[i], 
                    self.node[lineage[0]].state[1] + self.G0[idx, j],
                ))

        # Getting Current Length of Node
        len_ = len(node)

        for i in range(len(states)):
            self.node.append(nodeObj(
                states[i], 
                parent=lineage[0], 
                layer=len(lineage) + 1, 
                isTerminal= states[i][0] == self.finalBody,
                Δv = node[lineage[0]].Δv))

            # INCLUDING Δv CALCULATINGS AND TERMINAL STATES
            if self.node[-1].layer == 3:
                # self.debug()
                l0 = self.computeLambert(len(self.node)-1)
                Δv = self.computeΔv(len(self.node)-1, l0)
                self.node[-1].Δv -= Δv
                self.node[-1].Δv = max(0, self.node[-1].Δv)
                if self.node[-1].Δv == 0:
                    self.node[-1].isTerminal = True
            else:
                # print
                # print(self.node[-1].parent, len(self.node)-1)
                # print(self.spk.et2utc(self.node[lineage[0]].state[1], 'C', 14, 12), "→", self.spk.et2utc(self.node[-1].state[1], 'C', 14, 12))
                # self.debug()
                l0 = self.computeLambert(len(self.node)-1)
                l1 = self.computeLambert(lineage[0])
                Δv = self.computeΔv(len(self.node)-1, l0, l1)
                self.node[-1].Δv -= Δv
                self.node[-1].Δv = max(0, self.node[-1].Δv)
                if self.node[-1].Δv == 0:
                    self.node[-1].isTerminal = True

            if i == 0:
                id = len(node)

        # Assigning new children ID's to parent node
        self.node[lineage[0]].children = [i for i in range(len_, len(node))]

        return id

    def simulate(self, id):
        lineage = self.getLineage(id)
        parent = lineage[1]
        self.node[id].n += 1
        if self.node[id].Δv != 0:
            self.node[id].X = self.rnd.uniform(0, 1)
        else:
            self.node[id].X = 0    

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
        # print(Δv)
        return Δv

    def getStatesPath(self, id):
        states = []
        for id in self.getLineage(id)[:-1]:
            states.append(self.node[id].state)
        return states.reverse()

    def getResults(self):
        node = self.node
        attr = []
        head = ["ID", "Layer", "Lineage", "Δv Left"]
        for id in range(len(node)):
            if node[id].children == None and node[id].Δv != 0 and node[id].state[0] == self.finalBody:
                attr.append(
                    [id,
                    node[id].layer,
                    self.getLineage(id),
                    round(node[id].Δv, ndigits=5)
                    ]
                )
        print(self.tabulate.tabulate(attr, head))

    def plotPath(self, id):
        def axisEqual3D(ax):
            np = self.np
            extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
            sz = extents[:,1] - extents[:,0]
            centers = np.mean(extents, axis=1)
            maxsize = max(abs(sz))
            r = maxsize/2
            for ctr, dim in zip(centers, 'xyz'):
                getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

        plt = self.plt
        # STARTING PLOT WINDOW
        J2000_jd = 2451544.5
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.view_init(azim=0, elev=90)
        pk = self.pk
        lineage = self.getLineage(id)[:-1]
        for id in lineage:
            pl = pk.planet.jpl_lp(self.pkP[self.node[id].state[0]])
            pk.orbit_plots.plot_planet(
                pl, 
                t0 = pk.epoch(float(self.spk.et2utc(self.node[id].state[1], 'J', 10)[3:]) - J2000_jd), 
                axes = ax, 
                color=self.col[self.node[id].state[0]]
            )
        for id in lineage[:-1]:
            l = self.computeLambert(id)
            pk.orbit_plots.plot_lambert(l, axes = ax, color='c')
        axisEqual3D(ax)
        plt.show()

        """
        states = [('4', 737574186.1856555), ('2', 685636986.1856555), ('4', 668649666.1856555), ('3', 638971266.1856555)]
        states.reverse()
        """


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
                if i != 0:
                    state = (tree[i].state[0], self.spk.et2utc(tree[i].state[1], 'C', 14, 12))
                else:
                    state = []
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
                        state,
                        # type_,
                        tree[i].n,
                        tree[i].X,
                        tree[i].Δv,
                        tree[i].parent,
                        children,
                    ]
                )
            print(
                self.tabulate.tabulate(
                    cells,
                    ["id", "L", "Leaf", "State", "n", "X", "Δv", "Par", "Children"]
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
"""
from MCTS import *
mcts = MCTS('5', ["Apr 01, 2020", "Jun 01, 2020"], maxIters=5000, detail=16, oldTreeStyle=False)
mcts.getResults()
"""

if __name__ == "__main__":
    print("run")


