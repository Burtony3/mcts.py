"""
from MCTS import *
mcts = MCTS('7', ["Jan 01, 1997", "Jan 02, 1998"], maxIters=2000, detail=32, debug = True, C3max=70, flybyBodies = ['2', '3', '5'])
mcts = MCTS('5', ["May 01, 1989", "Jan 01, 1990"], maxIters=500, detail=8, debug = True, C3max=30)
mcts = MCTS('4', ["Jan 01, 2020", "Jan 01, 2022"], maxIters=3000, detail=16, debug=True)
mcts.getResults()
"""
from scipy.optimize import newton
class MCTS(object):
    def __init__(self): 

        # ------------------------- PACKAGE INSTALLATION ------------------------- #
        self.math     = __import__("math")
        self.np       = __import__("numpy")
        self.os       = __import__("os")
        self.pd       = __import__('pandas')
        self.plt      = __import__("matplotlib.pyplot")
        self.pk       = __import__("pykep")
        self.rnd      = __import__("random")
        self.spk      = __import__("spiceypy")
        self.sys      = __import__("sys")
        self.tabulate = __import__("tabulate")
        self.time     = __import__("time")
        # ------------------------------------------------------------------------ #

    def setup(
        self,
        finalBody,               # Target arrival body NAIF ID
        launchWindow,            # Beginning and End Date String (or single date)
        launchBody = '3',        # Launching body NAIF ID
        freeCapture = False,     # Ignores Δv to capture final body
        minimizeVinfCap = True,  # Whether to value high v∞ or low v∞ at final body 
        maxIters=10000,          # Maximum number of calculation iterations
        dvBudget=10,             # Maximum correction fuel budget
        detail=30,               # Number of indicies in launch & epoch nodes
        flybyBodies=['2', '3'],  # List of flyby body possibilities
        frame = 'ECLIPJ2000',    # Frame of reference for calculations
        abcorr = 'NONE',         # SPICE input
        cntrBody = '0',          # Central body for all state information
        C3max = 100,             # Maximum C3 capability of launch vehicle
        debug=False,             # Prints additional information
    ):

        # IMPORTING DSM .CSV
        self.dsmDvCsv = self.pd.read_csv("../matlab/dsm_integ/dsm_ega.csv", header=0)

        # SAVING SYSTEM CONSTRAINTS
        self.ΔvBudget     = dvBudget
        self.C3max        = C3max
        self.detail       = detail
        self.finalBody    = finalBody
        self.flybyBodies  = flybyBodies
        self.freeCapture  = freeCapture
        self.minVinfF     = minimizeVinfCap
        self.launchBody   = launchBody
        self.launchWindow = launchWindow
        self.maxIters     = maxIters

        # SAVING ENVIRONMENT PARAMENTERS
        self.frame = frame
        self.abcorr = abcorr
        self.cntrBody = cntrBody

        # SAVING CODE PARAMETERS
        self.debugBool = debug


        # FURNISHING KERNELS
        self.spk.furnsh("../data/spk/naif0009.tls")
        self.spk.furnsh("../data/spk/de438.bsp")
        self.spk.furnsh("../data/spk/gm_de431.tpc")
        self.spk.furnsh("../data/spk/pck00010.tpc")

        # INITIALIZING DICTIONARIES
        self.pkP = {'1': 'mercury', '2': 'venus',   '3': 'earth', '4': 'mars', '5': 'jupiter', '6': 'saturn',  '7': 'uranus',  '8': 'neptune'}
        self.tau = {'1': 87.97,     '2': 224.7,     '3': 365.25,  '4': 687,    '5': 4331,      '6': 10747,     '7': 30589,     '8': 59800}
        self.col = {'1': '#5a6165', '2': '#955c1c', '3': 'b',     '4': 'r',    '5': '#9e7a5b', '6': '#c9c0a9', '7': '#8eb0b8', '8': '#4d80d7'}

        # PARSING USER INPUTS
        self.P = flybyBodies + [finalBody]          # List of available planets for node states
        n = int(self.launchWindow[1][-4:]) - int(launchWindow[0][-4:])
        if len(launchWindow) == 2:                  # Checks for range or single launch window
            self.t0 = self.np.linspace(
                self.spk.utc2et(launchWindow[0]),
                self.spk.utc2et(launchWindow[1]),
                n*detail if n > 1 else detail       # Adds additional nodes if multiple years
            )
        else:
            self.t0 = self.spk.utc2et(launchWindow)

    def run(self):
        self.treeArt() # ASCII Graphic

        # SETTING STOP CONDITION (CTRL+C)
        import signal
        signal.signal(signal.SIGINT, self.debug)

        # STARTING LOG FILE
        if self.os.path.exists("log.txt"):
            self.os.remove("log.txt")

        # TAKING START TIME
        self.startTime = self.time.time()

        # CREATING TRACKING ELEMENTS
        self.runTime = 0
        self.numLam = 0

        # INITIALIZING TREE
        self.node = []                                                     # Allocating List
        self.node.append(nodeObj(None, parent = None, layer = 1))          # Creating Root Node
        self.node[0].children = []
        for i in range(len(self.t0)):                                      # Creating launch window children
            self.node.append(nodeObj((self.launchBody, self.t0[i]), layer = 2))
            self.node[0].children.append(len(self.node) - 1)

        # =============================================================== #
        # <><><><> RUNNING MONTE CARLO TREE SEARCH CREATION LOOP <><><><> #
        # =============================================================== #

        for itr in range(self.maxIters):
            if any([self.node[id].isTerminal for id in self.node[0].children]):
                print(itr)
                break
            # UPDATE ITERATION HEADER
            print("Current Iteration: {0}/{1}".format(itr, self.maxIters), end="\r")

            # PATH TO MOST VALUABLE LEAF NODE
            id = self.select() 

            # BREAK LOOP IF ALL NODES EXHAUSTED
            if id == -1:
                break

            # EXPANDS CURRENT LEAF IF VISITED PREVIOUSLY AND SELECTS CHILD
            if not self.node[id].children and (self.node[id].n != 0):
                id = self.expand(id)

            # ESIMATES FUTURE RETURNS FROM NODE
            X = self.simulate(id) 

            # UPDATES UP BRANCH
            self.backprop(id, X)

        self.runTime = self.time.time() - self.startTime

        # FINISHING PRINT
        if id != -1:
            print("Current Iteration: {0}/{0}".format(self.maxIters))
        else:
            print("Simulation Terminated All Possible Results Found")

        # DEBUGG PRINTING
        if self.debugBool:
            self.debug()

    ## ———————————————————————————————————————————————————————————————————————————— ##
    ## —————————————————————————————— SUB-ROUTINES ———————————————————————————————— ##
    ## ———————————————————————————————————————————————————————————————————————————— ##

    def select(self):
        def getReward(self, func, X, N, n):
            if func == 'ucb1':
                cp = 1 / self.math.sqrt(2.0)  # Cost Adjustment Parameter
                val = X + cp * self.math.sqrt(self.math.log1p(N) / n)
            elif func == 'epsG':
                print("FUNCTION NOT IMPLEMENTED")
                exit()
            else:
                print("NO FUNCTION DEFINED")
                exit()
            return val

        # CONSTANT INITIALIZATION
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
            for id_ in node[id].children:
                # GATHERING PROPERTIES
                X = node[id_].X                   # Getting child's future value
                N = node[id].n                    # Getting parent's number of visits
                n = node[id_].n                   # Getting child's number of visits
                isTerminal = node[id_].isTerminal # Checking if node is terminal/leaf

                if isTerminal:                           # Checking is node is terminal
                    ucb1.append(0)
                elif n == 0 and not isTerminal:          # If node is unvisited, immediately visit
                    id = id_
                    break
                else:                                    # Calculating Reward
                    ucb1.append(
                        getReward(self, 'ucb1', X, N, n)
                    )

            # CHECKING WHETHER UCB1 IS EMPTY
            if len(ucb1) is len_:
                # CHECKING IF ALL NODES ARE TERMINAL & RESTARTING SEARCH
                if all([val == 0 for val in ucb1]):
                    if id == 0: 
                        id = -1
                        self.logFile(id, "!!!All launch nodes considered terminal or fully visited")
                        return id
                    self.node[id].isTerminal = True
                    if self.node[id].layer < 5:
                        self.logFile(id, "!!Node {} in layer {} was considered terminal after creation".format(id, self.node[id].layer))
                    id = 0
                else:
                    indexChild = ucb1.index(max(ucb1)) # Finding index of child with max value
                    id = node[id].children[indexChild] # Finding its associated node ID

        return id

    def expand(self, id):
        # UPDATING VISITS
        self.node[id].n += 1

        # FIND IDs OF BRANCH
        lineage = self.getLineage(id)

        # NODE SHORT-HAND
        node = self.node

        # CREATING LIST OF STATES FOR NODE CREATION
        states = []
        P_ = self.flybyBodies + [self.finalBody]
        for i in range(len(P_)):
            et = self.setEpoch(self.node[id].state[1], self.node[id].state[0], P_[i])
            for j in range(len(et)):
                states.append((P_[i],et[j]))

        # FINDING CURRENT LENGTH OF ALL OF TREE
        len_ = len(node)


        for i in range(len(states)):
            self.node.append(nodeObj(
                states[i], 
                parent=id, 
                layer=len(lineage) + 1, 
                isTerminal= states[i][0] == self.finalBody,
                Δv = node[id].Δv))

            # INCLUDING Δv CALCULATINGS AND TERMINAL STATES
            if self.node[-1].layer == 3:                        # Launch Condition
                l0 = self.computeLambert(id = len(self.node)-1)
                dv = self.computeΔv(
                    id = len(self.node)-1, 
                    l0 = l0
                )
                vp = self.spk.spkezr(self.node[id].state[0], self.node[id].state[1], 'ECLIPJ2000', 'NONE', '0')[0][3:]
            else:                                               # Flyby Condition
                l0 = self.computeLambert(id = len(self.node)-1)
                l1 = self.computeLambert(id = lineage[0])
                dv, h = self.computeΔv(
                    id = len(self.node)-1, 
                    l0 = l0, 
                    l1 = l1
                )

            # ASSIGNING Δv TO NODE
            self.node[-1].Δv += dv
            self.node[-1].ΔvLeg = dv
            if self.node[-1].layer > 3:
                self.node[-1].h = h
            if self.node[-1].Δv > self.ΔvBudget or self.node[-1].h and self.node[-1].h < 0:
                self.node[-1].isTerminal = True
        
        # ASSIGN CREATED CHILDREN TO PARENT
        self.node[lineage[0]].children = [i for i in range(len_, len(node))]

        # FINDING HIGHEST VALUE CHILD TO SIMULATE
        id_ = []
        for i in self.node[lineage[0]].children:     # Collection all children not terminal
            if not self.node[i].isTerminal:
                id_.append(i)
        if len(id_) != 0:                            # Taking child with minimum Δv
            dv = [self.node[i].Δv for i in id_]
            id = id_[dv.index(min(dv))]
        if id == lineage[0]:                         # If all children terminal take first
            self.logFile(id, "!All children of node {} were considered terminal upon creation".format(id))
            id = self.node[lineage[0]].children[0]

        if all([self.node[i].isTerminal for i in self.node[lineage[0]].children]):
            self.node[lineage[0]].isTerminal = True
        
        return id

    def simulate(self, id):
        # FINDING IDs OF BRANCH
        lineage = self.getLineage(id)

        # UPDATING VISITS
        self.node[id].n += 1

        # COLLECTING LIST OF POSSIBLE PLANET STATES
        P_ = self.flybyBodies + [self.finalBody]

        # CREATING 2D LIST OF TEMPORARY STATES
        fakeStates = [(P, et) for P in P_ for et in self.setEpoch(self.node[id].state[1], self.node[id].state[0], P)]

        # INITIALIZING STATE EVALUATION CRITERIA
        X = []
        
        # RUNNING THROUGH LIST OF ALL PLAYABLE STATES
        for state in fakeStates:
            # SETTING INITIAL STATE
            state = [self.node[id].state, state]
            Δv = self.node[id].Δv
            l0 = self.computeLambert(id = id) if self.node[id].layer > 2 else None

            # RANDOM EXPLORATION UNTIL TERMINATION STATE
            i = 0
            while Δv < self.ΔvBudget or state[-1][0] != self.finalBody:
                i += 1
                if i == 10: print(state)
                # CALCULATING LAMBERT AND UPDATING Δv
                l1 = self.computeLambert(state0 = state[-2], state1 = state[-1])
                if len(state) == 2 and self.node[id].layer == 2:
                    Δv += self.computeΔv(id = id, l0 = l1, pState = "dum")
                else:
                    Δv += self.computeΔv(l0 = l0, l1 = l1, pState = state[-2])

                # APPENDING MATRIX AND RE-ASSIGNING VALUES
                l0 = l1
                P_ = self.flybyBodies + [self.finalBody]
                P_.remove(state[-1][0])
                i = self.rnd.randint(0, len(P_)-1)
                j = self.rnd.randint(0, self.detail-1)
                state.append((P_[i], self.setEpoch(state[-1][1], state[-1][0], P_[i])[j]))
                
            if Δv > self.ΔvBudget:
                bonus = 0.025
                X.append(bonus*(len(state)-2))
            else:
                X.append((self.ΔvBudget - Δv)/self.ΔvBudget)
                self.logFile(id, "Node {} child correctly simulated to final body with cost {}".format(id, X[-1]))

        X = sum(X)/len(X)

        return X

    def backprop(self, id, X):
        # FINDING IDs OF BRANCH
        lineage = self.getLineage(id)

        # AVERAGING IN NEW COST CALCULATED
        for id in lineage:
            self.node[id].X = (self.node[id].X*self.node[id].N + X)/(self.node[id].N + 1)
            self.node[id].N += 1

    ## -------------------------------------------------------------------------------- ##
    ## ------------------------------- HELPER FUNCTIONS ------------------------------- ##
    ## -------------------------------------------------------------------------------- ##

    def getLineage(self, id):
        """

        """
        lineage = []
        while id != 0:
            lineage.append(id)
            id = self.node[id].parent
        lineage.append(0)

        # Returns most recent node @ idx = 0 & oldest parent @ idx = -1
        return lineage

    def setEpoch(self, t0, p0, p1):
        # SETTING PARAMETERS FOR RANGES
        if p0 == p1: # If returning to same planet
            n = 1.05 # Return in ~2.1 periods
            m = 1.5  # Return in 3 periods
        elif int(p1) > 4:
            n = 0.05
            m = 0.5
        else:        # All other cases
            n = 0.1
            m = 1

        # SETTING EPOCH UPPER AND LOWER LIMITS
        et0 = n*(self.tau[p0] + self.tau[p1])*86400 + t0
        et1 = m*(self.tau[p0] + self.tau[p1] - 1)*86400 + t0

        # CREATING LINEAR RANGE
        et = self.np.linspace(et0, et1, self.detail)
        
        return et

    def queryDSM(self, dataframe, K, dTheta):
        pd = self.pd
        np = self.np
        # Current CSV is in 1 deg steps. Uncomment line below
        # when using this for final csv file (rounds to 0.5 deg)
        # dThetaLookup = round(dTheta*2)/2 

        # Placeholder Rounding dThetaLookup
        dThetaLookup = round(dTheta,0)

        rslt_df = dataframe[(dataframe['K'] == K) & (dataframe['Theta'] == dThetaLookup)]


        dsmdv  = np.float64(pd.to_numeric(rslt_df['DSMDV'])[0]).item()
        vinf1x = np.float64(pd.to_numeric(rslt_df['vinf1x'])[0]).item()
        vinf1y = np.float64(pd.to_numeric(rslt_df['vinf1y'])[0]).item()
        vinf   = np.array([vinf1x,vinf1y, 0.0])
        dt     = np.float64(pd.to_numeric(rslt_df['dt'])[0]).item()

        return dsmdv, vinf, dt

    def computeLambert(self, id = None, state0 = None, state1 = None):
        self.numLam += 1
        # IMPORTING TOOLBOXES
        spk = self.spk
        pk = self.pk

        if id:
            
            # GETTING PARENT
            parent = self.node[id].parent

            # GETTING STATES
            p0, et0 = self.node[parent].state
            p1, et1 = self.node[id].state
        else:
            p0, et0 = state0
            p1, et1 = state1

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

    def computeΔv(self, id = None, l0 = None, l1 = None, pState = None, oldDVMethod = False):
        # l₀ :: Lambert Arc from Parent → ID
        # l₁ ::                  Grandparent → Parent
        # IMPORTING REQUIRED PACKAGES
        spk = self.spk
        np = self.np
        math = self.math

        if id and not pState:
            parent = self.node[id].parent
            p = self.node[parent].state[0]
            et = self.node[parent].state[1]
        elif id and not l1 and type(pState) == str:
            p = self.node[id].state[0]
            et = self.node[id].state[1]
        else:
            p = pState[0]
            et = pState[1]
        vp = spk.spkezr(
            p,
            et,
            self.frame,
            self.abcorr,
            self.cntrBody
        )[0][3:]

        # CALCULATING Δv
        if not l1: # Returns Launch C3
            v0 = np.array(l0.get_v1()[0])/1000
            vinf = v0 - vp
            Δv = max(0, np.linalg.norm(vinf) - math.sqrt(self.C3max))
        else:
            # RETRIEVING TRAJECTORY VELOCITIES
            vi = np.array(l1.get_v2()[0])/1000 - vp
            vo = np.array(l0.get_v1()[0])/1000 - vp
            vih = vi + vp
            voh = vo + vp

            # GATHERING BODY DATA
            mu = spk.bodvrd(p + "99", "GM", 1)[1]
            R = spk.bodvrd(p + "99", "RADII", 3)[1][1]

            # SETTING NEWTON-RAPHSON INITIAL CONDITIONS
            aOutI = -mu / np.linalg.norm(vo)**2
            aOut = aOutI
            aIn = -mu/np.linalg.norm(vi)**2
            delta = math.acos(np.dot(voh, vih) / (np.linalg.norm(voh)*np.linalg.norm(vih)))

            # 0TH DERIVATIVE
            def f(eOut, aOut, aIn, delta):
                eOut = (aOut / aIn) * (eOut - 1) * math.sin(delta - math.asin(1 / eOut)) - 1
                return eOut

            # 1ST DERIVATIVE
            def f_prime(eOut, aOut, aIn, delta):
                tmp = [aOut / aIn, delta - math.asin(1 / eOut)]
                eOut = (tmp[0] * (eOut - 1) + 1) * (math.cos(tmp[1]) / (eOut**2 * math.sqrt(1 - eOut**-2))) + tmp[0] * math.sin(tmp[1])
                return eOut

            # CHECKING FOR FAIL CASE
            try:
                # ITERATING TO FIND OPTIMAL PERIAPSIS HEIGHT
                rp = aOutI * (1 - newton(f, 1.5, f_prime, (aOut, aIn, delta)))
                h = rp - R

                # FINDING Δv
                tmp = (2 * mu) / rp
                norm = np.linalg.norm
                Δv = abs(norm(math.sqrt(norm(vo)**2 + tmp)) - norm(math.sqrt(norm(vi)**2 + tmp)))
            except:
                Δv = 999.0
                h = -1

        # CHECKING IF ENTERING self.finalBody FOR CAPTURE Δv
        if id and self.node[id].state[0] == self.finalBody and not self.freeCapture:
            vp = self.spk.spkezr(
                self.finalBody,
                self.node[id].state[1],
                self.frame,
                self.abcorr,
                self.cntrBody
            )[0][3:] 
            vi = self.np.array(l0.get_v2()[0])/1000
            Δv += self.np.linalg.norm(vi - vp)

        # RETURN CASE
        if id and l1:       # Flyby
            return Δv, h
        else:               # Launch
            return Δv

    def getStatesPath(self, id):
        states = []
        for id in self.getLineage(id)[:-1]:
            states.append(self.node[id].state)
        return states.reverse()

    def getResults(self, minFlyby = 0, printR = True):
        node = self.node
        head = ["ID", "Flybys", "Lineage", "Δv Used", "V∞"]
        attr = []

        # FINDING NODES THAT MEET END CONDITION

        # PRUNING BASED ON RADIAL VELOCITY
        rm = []
        tmp = 3.0
        if int(self.finalBody) > 3:
            cond = lambda vr: vr > -tmp
        else:
            cond = lambda vr: vr < tmp

        for id in range(len(node)):
            if node[id].children == None and node[id].Δv < self.ΔvBudget and node[id].state[0] == self.finalBody and all([self.node[i].h > 0 for i in self.getLineage(id)[:-3]]):
                l = self.computeLambert(id)
                v = self.np.array(l.get_v2()[0])/1000
                vp = self.spk.spkezr(self.node[id].state[0], self.node[id].state[1], self.frame, self.abcorr, self.cntrBody)[0][3:]
                # print(self.np.array(l.get_v2()[0]), vp)
                vinf = self.np.linalg.norm(v - vp)

                # GETTING FINAL SPACECRAFT STATE
                r = self.spk.spkezr(self.node[id].state[0], self.node[id].state[1], self.frame, self.abcorr, self.cntrBody)[0][:3]
                vr = (self.np.dot(r, v)/self.np.linalg.norm(r))

                # CHECKING CONDITION
                if cond(vr):
                    attr.append(
                        [id,
                        node[id].layer-3,
                        self.getLineage(id),
                        round(node[id].Δv, 3),
                        round(vinf, 3)]
                    )

        # for id in [attr[i][0] for i in range(len(attr))]:
        #     # GETTING FINAL SPACECRAFT STATE
        #     l = self.computeLambert(id)
        #     r = self.spk.spkezr(self.node[id].state[0], self.node[id].state[1], self.frame, self.abcorr, self.cntrBody)[0][:3]
        #     v = self.np.array(l.get_v2()[0])/1000
        #     vr = (self.np.dot(r, v)/self.np.linalg.norm(r))

        #     # CHECKING CONDITION
        #     if not cond(vr):
        #         rm.append(id)

        # if rm:
        #     idL = [attr[i][0] for i in range(len(attr))]

        # SORTING BASED ON V∞ CONSTRAINT
        if self.minVinfF:
            attr = sorted(attr, key = lambda x: (1/x[3] if x[3] != 0 else float("inf"), x[1]), reverse=True)
        else:
            attr = sorted(attr, key = lambda x: (x[4], 1/x[3] if x[3] != 0 else float("inf"), x[1]), reverse=True)

        # PRINTING TO CONSOLE
        if printR: print(self.tabulate.tabulate(attr, head))

        # PULLING IDEAL NODES
        id_ = [attr[i][0] for i in range(len(attr))]
        return id_

    def plotPath(self, id, showPlanets = True, axes = None, lamC = "tab:blue", alpha = 1):
        """
        TODO: Add ability to create scrollable list
        """
        def axisEqual3D(ax):
            np = self.np
            extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
            sz = extents[:,1] - extents[:,0]
            centers = np.mean(extents, axis=1)
            maxsize = max(abs(sz))
            r = maxsize/2
            for ctr, dim in zip(centers, 'xyz'):
                getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

        plt = self.plt.pyplot

        pR = {'2': 0.7, '3': 1, '4':1.5, '5': 5.2, '6': 9.5, '7': 19.2, '8': 30}

        # STARTING PLOT WINDOW
        J2000_jd = 2451544.5
        if not axes:
            fig = plt.figure()
            ax = fig.gca(projection = '3d')
            ax.view_init(azim=0, elev=90)
        else:
            ax = axes
        pk = self.pk
        lineage = self.getLineage(id)[:-1]
        if showPlanets:
            for id in lineage:
                try:
                    pl = pk.planet.jpl_lp(self.pkP[self.node[id].state[0]])
                    ep = pk.epoch(float(self.spk.et2utc(self.node[id].state[1], 'J', 10)[3:]) - J2000_jd)
                    pk.orbit_plots.plot_planet(
                        pl, 
                        t0 = ep, 
                        axes = ax, 
                        color=self.col[self.node[id].state[0]],
                        units = 149597870700,
                        s = 0
                    )
                except:
                    theta = self.np.linspace(0, 2*self.np.pi, 100)
                    r = pR[self.node[id].state[0]]
                    z = self.np.linspace(0, 0, 100)
                    x = r*self.np.cos(theta)
                    y = r*self.np.sin(theta)
                    ax.plot3D(x, y, z, 'gray')

        for id in lineage[:-1]:
            l = self.computeLambert(id = id)
            pk.orbit_plots.plot_lambert(l, axes = ax, color=lamC, units = 149597870700, alpha = alpha)
        axisEqual3D(ax)
        if not axes:
            plt.show()

    def plotFamilies(self, idList, hideAxes = False, showTop = None):
        plt = self.plt.pyplot
        mpatches = self.plt.patches

        letterDict = {"2": "V", "3": "E", "4": "M", "5": "J", "6": "S", "7": "U", "8": "N"}

        # CREATING FIGURE WINDOW
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.view_init(azim=0, elev=90)

        cList = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan"]
        cList2 = ["b", "g", "r", "c", "m", "y", "k"]
        cList = ["tab:"+color for color in cList] + cList2
        cDict = {}

        # LOOPING THROUGH
        itr = 0
        if showTop and showTop < len(idList): idList = idList[:showTop]
        for id in idList:
            # GETTING NODE LINEAGE
            lineage = self.getLineage(id)[:-1]
            key = ''.join([self.node[i].state[0] for i in lineage[::-1]])

            # CHECKING FOR COLOR MATCHES KEY
            if key in cDict.keys():
                color = cDict[key]
            else:
                color = cList[itr]
                cDict[key] = color
                itr += 1

            # PLOTTING
            if id == idList[-1] or id == idList[0]:
                pBool = True
            else:
                pBool = False

            self.plotPath(id, showPlanets = pBool, axes = ax, lamC = color, alpha = 0.75)

        # ADDING FAMILIES TO LEGENDS
        patchList = []
        for key in cDict:
            dataKey = mpatches.Patch(color = cDict[key], label=''.join([letterDict[key[i]] for i in range(len(key))]))
            patchList.append(dataKey)

        plt.legend(handles = patchList)

        # TIGHTENING PLOT BOUNDARIES
        for spine in ax.spines.values():
            spine.set_visible(False)
        plt.tight_layout()

        # REMOVING AXES
        if hideAxes:
            plt.axis('off')
            plt.grid(b = None)

        plt.show()

    def export(self, filename = 'auto', exportFullTree = False):
        # IMPORTINT PACKAGES
        pd = self.pd
        spk = self.spk

        # FINDING RESUTLS
        id_ = self.getResults(printR = False)

        # STARTING XLSX
        if filename == 'auto':
            filename = "{0}_{1}_to_{2}_detail{3}_iters{4}.xlsx".format(
                self.pkP[self.finalBody], 
                self.launchWindow[0], 
                self.launchWindow[1], 
                self.detail, 
                self.maxIters
            )
        else:
            if filename[-5:] != '.xlsx':
                filename = filename + ".xlsx"
        writer = pd.ExcelWriter(filename, engine = 'xlsxwriter')

        # MCTS OUTPUT RESULTS
        out = {'id': [], 'Δv Total': [], 'C3': [], 'Lineage': [], 'Planets': [], 'Dates': [], 'Δv per Leg': [], 'Flyby Alt': [], 'Vinf': []}
        for id in id_:
            lineage = self.getLineage(id)[:-1] # Getting lineage minus root

            # STARTING ARRAYS
            P = []
            dates = []
            Δv = []
            h = []
            vinf = []

            # GETTING NODE STATES
            for i in lineage:
                P.append(self.node[i].state[0])
                dates.append(spk.et2utc(self.node[i].state[1], 'C', 14, 12))

            # GETTING LEG INFORMATION
            Δv = [round(self.node[i].ΔvLeg, 2) for i in lineage[:-1]]
            if not self.freeCapture:
                vp = self.spk.spkezr(P[0], self.node[id].state[1], self.frame, self.abcorr, self.cntrBody)[0][3:]
                l = self.computeLambert(id)
                vih = self.np.array(l.get_v2()[0])/1000
                vi = vih - vp
            else:
                vi = 0.0
            Δv[0] -= self.np.linalg.norm(vi)
            Δv.insert(0, round(self.np.linalg.norm(vi), 2))
            Δv[1] = round(Δv[1], 2)

            h = [round(float(self.node[i].h), 2) for i in lineage[:-2]]
                
            # CALCULATING C3
            vinf = self.node[lineage[-2]].Δv
            C3 = (self.math.sqrt(self.C3max) + vinf)**2


            l = self.computeLambert(id)
            v = self.np.array(l.get_v2()[0])/1000
            vp = self.spk.spkezr(self.node[id].state[0], self.node[id].state[1], self.frame, self.abcorr, self.cntrBody)[0][3:]
            # print(self.np.array(l.get_v2()[0]), vp)
            vinf = self.np.linalg.norm(v - vp)

            # APPENDING TO DICTIONARY
            out['id'].append(id)
            out['Δv Total'].append(round(self.node[id].Δv, 2))
            out['C3'].append(round(C3, 2))
            out['Lineage'].append(lineage[::-1])
            out['Planets'].append(P[::-1])
            out['Dates'].append(dates[::-1])
            out['Δv per Leg'].append(Δv[::-1])
            out['Flyby Alt'].append(h[::-1])
            out['Vinf'].append(vinf)

        # CONVERTING TO PANDAS & SAVING
        df = pd.DataFrame(out)
        df.to_excel(writer, sheet_name = 'Results')

        # EXPORTING FULL TREE
        if exportFullTree:
            out = {'id': [], 'visits': [], 'terminal': [], 'planet': [], 'date': [], 'dv': [], 'x': [], 'parent': [], 'children': [], 'vinf': []}
            for i in range(len(self.node)):
                # NODE INFORMATION
                out['id'].append(i)
                out['visits'].append(self.node[i].n)
                out['terminal'].append("True" if self.node[i].isTerminal else "False")

                # NODE STATE
                if self.node[i].state:
                    out['planet'].append(self.pkP[self.node[i].state[0]])
                    out['date'].append(spk.et2utc(self.node[i].state[1], 'C', 14, 12))
                else:
                    out['planet'].append("None")
                    out['date'].append("None")

                # NODE COST/REWARDS
                out['dv'].append(round(self.node[i].Δv, 2))
                out['x'].append(round(self.node[i].X, 4))

                # NODE "FAMILY"
                if self.node[i].parent:
                    out['parent'].append(self.node[i].parent)
                else:
                    out['parent'].append("None")
                out['children'].append(self.node[i].children)

                # NODE INBOUND V∞
                if self.node[i].parent and self.node[i].parent > self.detail:
                    l = self.computeLambert(i)
                    vp = self.spk.spkezr(self.node[i].state[0], self.node[i].state[1], self.frame, self.abcorr, self.cntrBody)[0][3:]
                    vinf = self.np.linalg.norm(self.np.array(l.get_v2()[0])/1000 - vp)
                    out['vinf'].append(vinf)
                else:
                    out['vinf'].append([])

            # EXPORTING TO SPREADSHEET
            df = pd.DataFrame(out)
            df.to_excel(writer, sheet_name = 'All Nodes')

        # GETTING RUN INFORMATION
        out = {
            'TREE SEARCH DIAGNOSITICS': [
            "RUN TIME:              {} SECONDS".format(round(self.runTime, 1)),
            "TOTAL OF LAMBERT RUNS: {}".format(self.numLam),
            "TOTAL NODES:           {}".format(len(self.node)),
            "ANGULAR DISPERSION:    {} DEGREES".format(360/self.detail)
        ]}
        df = pd.DataFrame(out)
        df.to_excel(writer, sheet_name = 'Run Information')

        writer.save()

    def logFile(self, id, statement):
        with open("log.txt", 'a') as f:
            f.write("{}\n".format(statement))


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

    def debug(self, *args):
        print(" ")
        print("### =================================================== ###")
        print("### ==================== DEBUGGING ==================== ###")
        print("### =================================================== ###")
        print(" ")
        time_ = self.time.time() - self.startTime
        cells = []
        tree = self.node
        if len(tree) > 500:
            r = input("PRINT ALL {0} NODES? (Y/N): ".format(len(tree)))
            r = r == "Y" or r == "y"
        else:
            r = True
        if r:
            print("↓TREE PRODUCED↓")
            length = 1001 if len(tree) > 1000 else len(tree)
            for i in range(length):
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
        print("RUN TIME:", round(time_, ndigits=5), " SECONDS")
        nLayer = []
        for i in range(len(tree)):
            if tree[i].layer > len(nLayer):
                nLayer.append(tree[i].n)
            else:
                nLayer[tree[i].layer - 1] += tree[i].n
        print("LAYER VISITS:", nLayer)
        numNodesInLayer = [0 for _ in range(len(nLayer))]
        for id in range(len(tree)): numNodesInLayer[tree[id].layer - 1] += 1
        print("NODES IN LAYERS:", numNodesInLayer)
        print("EQUIVALENT RESOLUTION PARAMETER:", 360/self.detail, "\bᵒ")
        print("NUMBER OF LAMBERT ARCS CALCULATED:", self.numLam)
        print(" ")

### ========================================================================== ###
### =============================== NODE CLASS =============================== ###
### ========================================================================== ###


class nodeObj:
    def __init__(self, state, parent=0, Δv=0, vInfIn=float("NaN"), layer=None, isTerminal=False):
        self.n = 0  # Visits in UCB1 Function
        self.N = 0
        # Estimated Reward in UCB1 Function (if -1 = dead node?)
        self.X = 0
        self.Δv = Δv  # Δv used up to that point
        self.ΔvLeg = 0
        self.h = []
        # Establishes Parent above (Single Scalar of IDs)
        self.parent = parent
        self.children = None  # Children Below (Array of IDs)
        self.state = state  # Either NAIF ID String or Epoch
        self.layer = layer
        self.isTerminal = isTerminal