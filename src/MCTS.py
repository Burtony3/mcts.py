from scipy.optimize import newton
class MCTS(object):
    def __init__(self): 

        # ------------------------- PACKAGE INSTALLATION ------------------------- #
        self.click    = __import__("click")
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

    def loadKernels(self):
        self.spk.furnsh(self.filepath + "/../data/spk/naif0009.tls")
        self.spk.furnsh(self.filepath + "/../data/spk/de438.bsp")
        self.spk.furnsh(self.filepath + "/../data/spk/gm_de431.tpc")
        self.spk.furnsh(self.filepath + "/../data/spk/pck00010.tpc")

    def setup(
        self,
        finalBody,               # Target arrival body NAIF ID
        launchWindow,            # Beginning and End Date String (or single date)
        launchBody = '3',        # Launching body NAIF ID
        freeCapture = False,     # Ignores Δv to capture final body
        minimizeVinfCap = True,  # Whether to value high v∞ or low v∞ at final body 
        maxIters=10000,          # Maximum number of calculation iterations
        dvBudget=10,             # Maximum correction fuel budget
        detail=16,               # Number of indicies in launch & epoch nodes
        flybyBodies=['2', '3'],  # List of flyby body possibilities
        enableDSM = True,        # Allow for Earth Targeting DSMs
        frame = 'ECLIPJ2000',    # Frame of reference for calculations
        abcorr = 'NONE',         # SPICE input
        cntrBody = '0',          # Central body for all state information
        C3max = 100,             # Maximum C3 capability of launch vehicle
        debug=False,             # Prints additional information
    ):

        # DEFINING CONSTANTS
        self.AU2m = 149597870700
        self.AU2km = self.AU2m/1000

        # SAVING SYSTEM CONSTRAINTS
        self.ΔvBudget     = dvBudget
        self.C3max        = C3max
        self.detail       = detail
        self.enableDSM    = enableDSM
        self.finalBody    = finalBody
        self.flybyBodies  = flybyBodies
        self.freeCapture  = freeCapture
        self.minVinfF     = minimizeVinfCap
        self.launchBody   = launchBody
        self.launchWindow = launchWindow
        self.maxIters     = maxIters
        self.filepath     = self.os.getcwd()
        
        # IMPORTING DSM .CSV
        self.dsmDvCsv = self.pd.read_csv(self.filepath + "/../matlab/dsm_integ/dsm_ega.csv", header=0)

        # SAVING ENVIRONMENT PARAMENTERS
        self.frame = frame
        self.abcorr = abcorr
        self.cntrBody = cntrBody

        # SAVING CODE PARAMETERS
        self.debugBool = debug

        # FURNISHING KERNELS
        self.loadKernels()

        # INITIALIZING DICTIONARIES
        self.pkP = {'1': 'mercury', '2': 'venus',   '3': 'earth', '4': 'mars', '5': 'jupiter', '6': 'saturn',  '7': 'uranus',  '8': 'neptune'}
        self.tau = {'1': 87.97,     '2': 224.7,     '3': 365.25,  '4': 687,    '5': 4331,      '6': 10747,     '7': 30589,     '8': 59800}
        self.col = {'1': '#5a6165', '2': '#955c1c', '3': 'b',     '4': 'r',    '5': '#9e7a5b', '6': '#c9c0a9', '7': '#8eb0b8', '8': '#4d80d7'}
        self.letterDict = {"2": "V", "3": "E", "4": "M", "5": "J", "6": "S", "7": "U", "8": "N"}

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

    def run(self, iters = None, startNode = 0):
        self.treeArt() # ASCII Graphic

        # SETTING STOP CONDITION (CTRL+C)
        self.breakLoop = False
        import signal
        signal.signal(signal.SIGINT, self.debug)
        signal.signal(signal.SIGTSTP, self.exit)

        # STARTING LOG FILE
        if not iters:
            self.startNode = 0
            maxIters = self.maxIters
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
        else:
            self.startNode = startNode
            self.maxIters += iters
            maxIters = iters

        # =============================================================== #
        # <><><><> RUNNING MONTE CARLO TREE SEARCH CREATION LOOP <><><><> #
        # =============================================================== #

        with self.click.progressbar(range(maxIters), fill_char = '█', empty_char = ' ', label='Building Trajectories') as bar:
            for itr in bar:

                # CHECKING FOR BREAK CONDITION
                if all([self.node[id].isTerminal for id in self.node[self.startNode].children]):
                    print(itr)
                    id = -1
                    break

                # PATH TO MOST VALUABLE LEAF NODE
                id = self.select() 

                # BREAK LOOP IF ALL NODES EXHAUSTED
                if id == -1:
                    break

                # EXPANDS CURRENT LEAF IF VISITED PREVIOUSLY AND SELECTS CHILD
                if not self.node[id].children and self.node[id].n != 0:
                    id = self.expand(id)

                # ESIMATES FUTURE RETURNS FROM NODE
                X = self.simulate(id) 

                # UPDATES UP BRANCH
                self.backprop(id, X)

                if self.breakLoop:
                    break

            # RECORDING RUN TIME
            self.runTime = self.time.time() - self.startTime

    ## ———————————————————————————————————————————————————————————————————————————— ##
    ## —————————————————————————————— SUB-ROUTINES ———————————————————————————————— ##
    ## ———————————————————————————————————————————————————————————————————————————— ##

    def select(self):
        # REWARD FUNCTION
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
        id = self.startNode           # Starting at top of tree
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
                # CHECKING IF ALL NODES ARE TERMINAL & RESTARTING SEARCHΔ
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

        # DICT FOR DSM θ RANGES
        kDict = {
            '2-': [ -7,   -75  ],
            '2+': [  7,    63  ],
            '3-': [ -7,   -75  ],
            '3+': [ 13,    42.5], #← Holes in CSV
            '4-': [-31.5, -75  ],
            '4+': [ 26,    42.5]
        }

        # INITIALIZING STATE PAIRS
        states = []
        P_ = self.flybyBodies + [self.finalBody]

        # PRUNING CERTAIN SEQUENCES
        # REMOVING EARTH FLYBY's WITH DSMs
        if len(self.node[id].state[0]) != 3 and self.enableDSM and self.node[id].state[0] in P_: 
            P_.remove(self.node[id].state[0])

        # REMOVING SEQUENCES OF MORE THAN 2 OF THE SAME PLANET IN A ROW
        if self.node[id].state[0] in P_ and self.node[id].layer > 2 and self.node[id].state[0] == self.node[lineage[1]].state[0]:
            P_.remove(self.node[id].state[0])

        # ADDING DSMs
        if self.node[id].layer == 2 and self.enableDSM:
            P_.append('32-')
            P_.append('32+')
            P_.append('33-')
            P_.append('33+')
            P_.append('34-')
            P_.append('34+')


        for P__ in P_:
            if "+" in P__ or "-" in P__:                                    # DSM Case
                et = [0 for _ in range(int(self.detail*0.5))]
                θ = []
                lim = kDict[P__[1:]]
                for idx in range(int(self.detail*0.5)):
                    tmp = lim[0] + (lim[1] - lim[0])*((idx)/(int(self.detail*0.5)))
                    θ.append(round(tmp*2)/2)
            else:                                                    # Normal Case
                et = self.setEpoch(self.node[id].state[1], self.node[id].state[0], P__)
            for i in range(len(et)):
                if "+" in P__ or "-" in P__:                   # DSM Case
                    states.append((P__, et[i], θ[i]))
                else:                                   # Normal Case
                    states.append((P__, et[i]))

        # FINDING CURRENT LENGTH OF ALL OF TREE
        len_ = len(node)

        # CREATING NODES BASED OFF STATE PAIRS
        for i in range(len(states)):
            self.node.append(nodeObj(
                states[i], 
                parent=id, 
                layer=len(lineage) + 1, 
                isTerminal= states[i][0] == self.finalBody,
                Δv = node[id].Δv))

            # INCLUDING Δv CALCULATINGS AND TERMINAL STATES
            if self.node[-1].layer == 3:                            # Launch Condition
                delta = None
                if "+" not in self.node[-1].state[0] and "-" not in self.node[-1].state[0]:   # Normal Case
                    l0 = self.computeLambert(id = len(self.node)-1)
                    dv = self.computeΔv(
                        id = len(self.node)-1, 
                        l0 = l0
                    )
                else:                                               # DSM Condition
                    dv, _, dt, vinfL = self.queryDSM(self.node[-1].state[0][1], self.node[-1].state[2])
                    self.node[-1].state = (self.node[-1].state[0], self.node[id].state[1] + dt*86400, self.node[-1].state[2])
                    dv += max(0, vinfL - self.math.sqrt(self.C3max))
            else:                                                   # Flyby Condition
                l0 = self.computeLambert(id = len(self.node)-1)
                l1 = self.computeLambert(id = lineage[0])
                dv, h, delta = self.computeΔv(
                    id = len(self.node)-1, 
                    l0 = l0, 
                    l1 = l1
                )

            # ASSIGNING Δv TO NODE
            self.node[-1].Δv += dv
            self.node[-1].ΔvLeg = dv

            # SAVING FLYBY VARIABLES IF APPLICABLE
            if delta:
                self.node[-1].delta = self.np.degrees(delta)
            if self.node[-1].layer > 3:
                self.node[-1].h = h

            # UPDATING NODE TERMINALITY
            if self.node[-1].Δv > self.ΔvBudget or self.node[-1].h and self.node[-1].h < 150:
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
            # print(self.node[id].state)
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

        X = sum(X)/len(X) * 10 # Amplifying scores to increase gap

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
        lineage = []
        while id != 0:
            lineage.append(id)
            id = self.node[id].parent
        lineage.append(0)

        # Returns most recent node @ idx = 0 & oldest parent @ idx = -1
        return lineage

    def setEpoch(self, t0, p0, p1):
        # SETTING PARAMETERS FOR RANGES
        if p0 == p1:      # If returning to same planet
            det = self.detail // 3
            et = self.np.array([])
            for K in range(1, 4):
                et0 = 0.9  * K * self.tau[p0[0]]*86400 + t0
                et1 = (K * self.tau[p0[0]] - 3)*86400 + t0
                et_ = self.np.linspace(et0, et1, det)
                for val in et_: et = self.np.append(et, val)
        
        else:
            if int(p1) > 4:   # Outer Planets Case
                n = 0.05
                m = 0.25
            else:             # All other cases
                n = 0.1
                m = 1

            # SETTING EPOCH UPPER AND LOWER LIMITS
            et0 = n*(self.tau[p0[0]] + self.tau[p1[0]])*86400 + t0
            et1 = m*(self.tau[p0[0]] + self.tau[p1[0]] - 3)*86400 + t0

            # CREATING LINEAR RANGE
            et = self.np.linspace(et0, et1, self.detail)
        
        return et

    def queryDSM(self, K, dThetaLookup):
        pd = self.pd
        np = self.np

        # PULLING DATAFRAME
        dataframe = self.dsmDvCsv

        # FINDING CORRECT ROW ID
        rslt_df = dataframe[(dataframe['K'] == int(K)) & (dataframe['Theta'] == dThetaLookup)]

        # QUERYING INFORMATION
        dsmdv  = pd.to_numeric(rslt_df['DSMDV']).item()
        vinfL  = pd.to_numeric(rslt_df['vinfLaunch']).item()
        vinf1x = pd.to_numeric(rslt_df['vinf1x']).item()
        vinf1y = pd.to_numeric(rslt_df['vinf1y']).item()
        vinf   = np.array([vinf1x, vinf1y])
        dt     = pd.to_numeric(rslt_df['dt']).item()

        return dsmdv, vinf, dt, vinfL

    def computeLambert(self, id = None, state0 = None, state1 = None):
        self.numLam += 1
        # IMPORTING TOOLBOXES
        spk = self.spk
        pk = self.pk

        if id:                                  # Normal Case
            # GETTING PARENT
            parent = self.node[id].parent

            # GETTING STATES
            p0  = self.node[parent].state[0][0]
            et0 = self.node[parent].state[1]
            p1  = self.node[id].state[0][0]
            et1 = self.node[id].state[1]
        else:                                   # Random Exploration Case
            p0  = state0[0][0]
            et0 = state0[1]
            p1  = state1[0][0]
            et1 = state1[1]

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

        if id and not pState:                       # Flyby Condition
            parent = self.node[id].parent
            p = self.node[parent].state[0][0]
            et = self.node[parent].state[1]
            # if len(p) > 
        elif id and not l1 and type(pState) == str: # Launch C3 Condition
            p = self.node[id].state[0][0]
            et = self.node[id].state[1]
        else:                                       # Random Exploration Condition
            p = pState[0][0]
            et = pState[1]

        # FINDING FLYBY/LAUNCH PLANET VELOCITY
        vp = spk.spkezr(
            p,
            et,
            self.frame,
            self.abcorr,
            self.cntrBody
        )[0][3:]

        # CALCULATING Δv
        # REQUIRED ADDITIONAL C3 Δv
        if not l1:
            v0 = np.array(l0.get_v1()[0])/1000
            vinf = v0 - vp
            Δv = max(0, np.linalg.norm(vinf) - math.sqrt(self.C3max))

        # FLYBY Δv
        else:
            # FINDING V∞i
            # FROM DSM CASE
            if id and not pState and len(self.node[parent].state) == 3:
                # FINDING INITIAL DSM LAUNCH PLANET
                gparent = self.node[parent].parent

                # FINDING ENTRY VELOCITY STATE INTO FLYBY PLANET
                _, vi, _, _ = self.queryDSM(self.node[parent].state[0][1], self.node[parent].state[2])

                # ROTATE VECTOR TO ALIGN WITH J2000 FRAME FROM ROTATEING FRAME
                pE = self.spk.spkezr(self.node[gparent].state[0], self.node[gparent].state[1], self.frame, self.abcorr, self.cntrBody)[0][:3]
                theta = self.getJ2000Ang(pE)
                theta = 2*np.pi - (theta + self.np.radians(90))
                c, s = (self.np.cos(theta), self.np.sin(theta))
                R = np.array([[c, -s], [s, c]])
                vi = self.np.dot(vi, R)
                vi = self.np.append(vi, 0)

            # FROM LAMBERT ARC CASE
            else:
                vi = np.array(l1.get_v2()[0])/1000 - vp

            # OUTBOUND VELOCITY
            vo = np.array(l0.get_v1()[0])/1000 - vp

            # GATHERING BODY DATA
            mu = spk.bodvrd(p + "99", "GM", 1)[1]
            R = spk.bodvrd(p + "99", "RADII", 3)[1][1]

            # SETTING NEWTON-RAPHSON INITIAL CONDITIONS
            aOutI = -mu / np.linalg.norm(vo)**2
            aOut = aOutI
            aIn = -mu/np.linalg.norm(vi)**2
            delta = self.math.acos(self.np.dot(vi, vo)/(self.np.linalg.norm(vi)*self.np.linalg.norm(vo)))

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
            return Δv, h, delta
        else:               # Launch
            return Δv

    def getStatesPath(self, id):
        states = []
        for id in self.getLineage(id)[:-1]:
            states.append(self.node[id].state)
        return states.reverse()

    def getResults(self, minFlyby = 0, minVr = 3, printR = False):
        node = self.node
        head = ["ID", "Flybys", "Lineage", "Δv Used", "V∞"]
        attr = []

        # FINDING NODES THAT MEET END CONDITION

        # PRUNING BASED ON RADIAL VELOCITY
        if int(self.finalBody) > 3:
            cond = lambda vr: vr > minVr
        else:
            cond = lambda vr: vr < minVr

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

    def axisEqual3D(self, ax):
        np = self.np
        extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
        sz = extents[:,1] - extents[:,0]
        centers = np.mean(extents, axis=1)
        maxsize = max(abs(sz))
        r = maxsize/2
        for ctr, dim in zip(centers, 'xyz'):
            getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

    def getJ2000Ang(self, r):
        theta = self.math.acos(r[0]/self.np.linalg.norm(r))
        if r[1] < 0:
            theta = -theta
        return theta

    def plotPath(self, id, showPlanets = True, axes = None, lamC = "tab:blue", alpha = 1, cntrS = 50):
        # IMPORTING PACKAGES
        pk = self.pk
        plt = self.plt.pyplot

        # STARTING PLOT WINDOW
        if not axes:
            fig = plt.figure()
            ax = fig.gca(projection = '3d')
            ax.view_init(azim=0, elev=90)
        else:
            ax = axes

        # RETRIEVING LINEAGE
        lineage = self.getLineage(id)[:-1]

        # PLOTTING LAMBERT ARCS
        for id in lineage[:-1]:
            # DSM FLYBY CASE (NOT ACTUAL ORBIT)
            if len(self.node[id].state) == 3:
                # RETRIEVING PARENT AND STATE
                parent = self.node[id].parent
                X = self.spk.spkezr(self.node[parent].state[0][0], self.node[parent].state[1], self.frame, self.abcorr, self.cntrBody)[0]

                # RETRIVING DSM INFORMATION
                K = self.node[id].state[0][1]
                _, _, dt, vinfL = self.queryDSM(K, self.node[id].state[2])

                # CALCULATING DSM APEHELION
                rp1 = X[:3]                         # Perihelion position
                theta = self.getJ2000Ang(rp1)       # Angle of perihelion wrt J2000
                vp = X[3:]                          # Perihleion velocity
                r = self.np.linalg.norm(rp1)        # Normal Value
                v = self.np.linalg.norm(vp) + vinfL # Normal Value
                mu = self.spk.bodvrd('sun', 'GM', 1)[1].item() # μ of Sun
                a = -mu/(v**2 - (2*mu)/r)           # Semi-Major Axis
                tof = self.math.pi*self.math.sqrt(a**3 / mu) # ½ period time of flight
                ra = 2*a - r                        # Apehelion normal value           
                phi = theta + self.math.pi + self.np.radians(1) # Angle of apehelion wrt J2000 (+ small bias)
                ra = [ra*self.math.cos(phi), ra*self.math.sin(phi), 0]  # Apehelion vector

                # MARKING DSM MANEUVER POINT
                ax.scatter(ra[0]/self.AU2km, ra[1]/self.AU2km, ra[2]/self.AU2km, color = lamC, marker='x')

                # CREATING LAMBERT ARC INPUTS
                ra = self.np.array(ra)
                rp2 = self.spk.spkezr(self.node[id].state[0][0], self.node[id].state[1], self.frame, self.abcorr, self.cntrBody)[0][:3]
                
                # LAMBERT ARC FROM LAUNCH TO DSM POINT
                l1 = self.pk.lambert_problem(
                    r1 = rp1*1000,
                    r2 = ra*1000,
                    tof = tof,
                    mu = self.pk.MU_SUN
                )
                pk.orbit_plots.plot_lambert(l1, axes = ax, color=lamC, units = self.AU2m, alpha = alpha)

                # LAMBERT ARC FROM DSM POINT TO FLYBY
                l2 = self.pk.lambert_problem(
                    r1 = ra*1000,
                    r2 = rp2*1000,
                    tof = (dt*86400 - tof),
                    mu = self.pk.MU_SUN
                )
                pk.orbit_plots.plot_lambert(l2, axes = ax, color=lamC, units = self.AU2m, alpha = alpha)

            # NORMAL LAMBERT ARC CASE
            else:
                l = self.computeLambert(id = id)
                pk.orbit_plots.plot_lambert(l, axes = ax, color=lamC, units = self.AU2m, alpha = alpha)

        # PLOTTING PLANETS AVAILABLE FOR FLYBY AND ARRIVAL
        if showPlanets:            
            # FINDING INITIAL EPOCH
            et0 = self.node[id].state[1]
            if self.enableDSM: self.P.append("3") # Adding Earth for DSM Case
            for P in self.P:
                # FINDING SECONDS IN 1 ORBIT PERIOD
                tau = self.tau[P]*86400
                et1 = et0 + tau
                et = self.np.linspace(et0, et1, 250)

                # COLLECTING DATA POINTS
                x = []
                y = []
                z = []
                for et_ in et:
                    r = self.spk.spkezr(P, et_, self.frame, self.abcorr, self.cntrBody)[0][:3]/self.AU2km
                    x.append(r[0])
                    y.append(r[1])
                    z.append(r[2])
                ax.plot3D(x, y, z, 'black', linewidth=0.5)

            # PLOTTING CENTER BODY SURFACE
            ax.scatter3D(0, 0, 0, s = cntrS, c = "#FFCC33")

        if not axes:
            self.axisEqual3D(ax)
            ax.set_box_aspect((1, 1, 1))
            plt.show()

    def plotFamilies(self, idList, hideAxes = False, showTop = None):
        print("CHANGED")
        # IMPORTING PACKAGES
        plt = self.plt.pyplot
        mpatches = self.plt.patches


        # CREATING FIGURE WINDOW
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.view_init(azim=0, elev=90)

        # CREATING COLOR LISTS
        cList = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan"]
        cList2 = ["b", "g", "r", "c", "m", "y", "k"]
        cList = ["tab:"+color for color in cList] + cList2
        tmpDict = {"2": "a", "3": "b", "4": "c"}

        # INITIALIZING DICTIONARY RELATING SEQUENCE WITH PLOT COLOR
        cDict = {}  

        # LOOPING THROUGH
        itr = 0
        if showTop and showTop < len(idList): idList = idList[:showTop]
        for id in idList:
            # GETTING NODE LINEAGE
            lineage = self.getLineage(id)[:-1]
            key = []

            # CREATING SEQUENCE STRINGS (IE: EVEEJ)
            for i in lineage[::-1]:
                if "+" in self.node[i].state[0] or "-" in self.node[i].state[0]:
                    K = self.node[i].state[0][1]
                    p = self.node[i].state[0][2]
                    key.append("({}:1{}Δv)E".format(tmpDict[K], p))
                else:
                    key.append(self.node[i].state[0])
            key = ''.join(key)

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
            lab = []
            for i in range(len(key)):
                    if key[i] in tmpDict.values():
                        lab.append([key_ for key_ in tmpDict.items() if key_[1] == key[i]][0][0])
                    else:
                        lab.append(self.letterDict[key[i]] if key[i] in self.letterDict else key[i])
            lab = ''.join(lab)
            dataKey = mpatches.Patch(color = cDict[key], label=lab)
            patchList.append(dataKey)

        plt.legend(handles = patchList)

        # REMOVING AXES
        if hideAxes:
            plt.axis('off')
            plt.grid(b = None)
        else:
            ax.set_xlabel("X Position (AU)")
            ax.set_ylabel("Y Position (AU)")
            ax.set_zlabel("Z Position (AU)")

        # TIGHTENING PLOT BOUNDARIES
        # for spine in ax.spines.values():
            # spine.set_visible(False)
        # plt.tight_layout()

        ax.set_title("Trajectories to {}\n{} to {}\nΔv Max: {} km/s".format(self.pkP[self.finalBody].capitalize(), self.launchWindow[0], self.launchWindow[1], self.ΔvBudget))

        self.axisEqual3D(ax)
        ax.set_box_aspect((1, 1, 1))

        plt.show()

    def queryBranchAttributes(self, id):
        # GETTING RID OF NUMPY WARNING
        import warnings
        warnings.simplefilter(action='ignore', category=FutureWarning)

        # GETTING NODE LINEAGE
        lineage = self.getLineage(id)[-2::-1] # Excluding Root

        # INITIALIZING DICTIONARY
        attr = {}

        # RETRIEVING BRANCH SPECIFIC INFORMATION
        tof   = (self.node[lineage[-1]].state[1] - self.node[lineage[0]].state[1])/86400
        dvTot = self.node[id].Δv

        # GETTING LAUNCH C3
        if len(self.node[lineage[1]].state) == 3:
            _, _, _, vinfL = self.queryDSM(self.node[lineage[1]].state[0][-2], self.node[lineage[1]].state[2])
            C3 = vinfL**2
        else:
            l  = self.computeLambert(lineage[1])
            v  = self.np.array(l.get_v1()[0])/1000
            vp = self.spk.spkezr(self.node[lineage[0]].state[0], self.node[lineage[0]].state[1], self.frame, self.abcorr, self.cntrBody)[0][3:]
            C3 = self.np.linalg.norm(v - vp)**2

        # GETTING ARRIVAL V∞
        l     = self.computeLambert(id)
        v     = self.np.array(l.get_v2()[0])/1000
        if len(self.node[id].state) == 3:
            vp    = self.spk.spkezr(self.node[id].state[0][:-2], self.node[id].state[1][:-2], self.frame, self.abcorr, self.cntrBody)[0][3:]
        else:
            vp    = self.spk.spkezr(self.node[id].state[0], self.node[id].state[1], self.frame, self.abcorr, self.cntrBody)[0][3:]
        vInfF = self.np.linalg.norm(v - vp)

        # RETRIEVING NODE SPECIFIC INFORMATION
        P      = [self.node[id].state[0][:-2] if len(self.node[id].state) == 3 else self.node[id].state[0] for id in lineage]
        utc    = [self.spk.et2utc(self.node[id].state[1], 'C', 14, 12) for id in lineage]
        theta  = [self.node[id].state[2] if len(self.node[id].state) == 3 else "None" for id in lineage]
        K      = [self.node[id].state[0][-2] if len(self.node[id].state) == 3 else "None" for id in lineage]
        n      = [self.node[id].n for id in lineage]
        dateL  = self.spk.utc2et(utc[0])

        # RETRIEVING LEG SPECIFIC INFORMATION
        lineage = lineage[1:]  # Looking at every Leg
        dvLeg   = [self.node[id].ΔvLeg for id in lineage]
        if not self.freeCapture: dvLeg[-1] = dvLeg[-1] - vInfF
        if not self.freeCapture: dvLeg.append(vInfF)
        h       = [self.node[id].h.item() if self.node[id].h else "None" for id in lineage]
        delta   = [self.node[id].delta if self.node[id].delta else "None" for id in lineage]

        # GETTING LEG V∞
        vInfOVec = []
        vInfIVec = []
        for id in lineage:
            # FINDING LAMBERT
            l = self.computeLambert(id)

            # FINDING OUTPUT V∞ OF LEG
            vp = self.spk.spkezr(
                self.node[self.node[id].parent].state[0][:-2] if len(self.node[self.node[id].parent].state) == 3 else self.node[self.node[id].parent].state[0], 
                self.node[self.node[id].parent].state[1], 
                self.frame, self.abcorr, self.cntrBody
            )[0][3:]
            vInfOVec.append(self.np.array(l.get_v1()[0])/1000 - vp)

            # FINDING INBOUND V∞ OF LEG
            vp = self.spk.spkezr(
                self.node[id].state[0][:-2] if len(self.node[id].state) == 3 else self.node[id].state[0], 
                self.node[id].state[1], 
                self.frame, self.abcorr, self.cntrBody
            )[0][3:]
            vInfIVec.append(self.np.array(l.get_v2()[0])/1000 - vp)

        # ADJUSTING LIST TO SHOW FOR EACH PLANET
        vInfOVec.append("None")
        vInfIVec.insert(0, "None")

        # FINDING NORMAL VALUES
        vInfO = [self.np.linalg.norm(val) if val != "None" else "None" for val in vInfOVec]
        vInfI = [self.np.linalg.norm(val) if val != "None" else "None" for val in vInfIVec]
            

        # CREATING SEQUENCE NAME
        seq = []
        for id in range(len(P)):
            if K[id] != "None":
                seq.append("({}:1{}){}".format(K[id], "-" if str(self.np.sign(theta[id]))[0] == "-" else "+", self.letterDict[P[id]]))
            else:
                seq.append(self.letterDict[P[id]])
        seq = ''.join(seq)

        # SAVING TO DICTIONARY
        varSkip = ["lineage", "self", "attr", "id", "v", "vp", "l", "vinfL", "varSkip", "warnings", "_"]
        for name in dir():
            if not name.startswith("__") and name not in varSkip:
                attr[name] = eval(name)

        return attr


    def plotAttributes(self, idList, xAttr, yAttr, cAttr, cmap = 'turbo'):
        # IMPORTING TOOLBOXES
        plt = self.plt.pyplot

        # INITIALIZING FIGURE
        fig = plt.figure()
        ax = fig.gca()

        # SETTING LABELS & TITLES
        ax.set_xlabel(xAttr)
        ax.set_ylabel(yAttr)
        ax.set_title("To {}\n{} to {}\n{}k Iterations".format(self.pkP[self.finalBody].capitalize(), self.launchWindow[0], self.launchWindow[1], self.maxIters/1000))

        # INITIALIZING PLOTTED VALUES
        x = []
        y = []
        c = []
        seq = []
        Δv = []

        # RETRIEVING PLOTTED VALUES
        for id in idList:
            attr = self.queryBranchAttributes(id)
            x.append(attr[xAttr])
            y.append(attr[yAttr])
            c.append(attr[cAttr])
            seq.append(attr['seq'])
            Δv.append(attr['dvTot'])

        # PLOTTING
        ax.grid('on', zorder = 3)
        sc = ax.scatter(x, y, c = c, s=7.5, zorder = 3, cmap=cmap)

        # ADDING COLORBAR
        cbar = fig.colorbar(sc)
        cbar.ax.set_ylabel(cAttr)

        # LAUNCH DATE EDGE CASE
        if xAttr == "dateL":
            vals, _ = plt.xticks()
            labels = [self.spk.et2utc(val, 'C', 14, 12) for val in vals]
            plt.xticks(vals, labels, rotation = 30, ha = 'right', size='small')

        # CREATING ANNOTATIONS
        annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)

        def update_annot(ind):
            idx = ind["ind"][0]
            pos = sc.get_offsets()[idx]
            annot.xy = pos
            """
            l, r = plt.xlim()
            if pos[0] > l + 0.5*(r - l):
                annot.xytext = (-annot.xytext[0], annot.xytext[1])
            b, t = plt.ylim()
            if pos[1] > b + 0.5*(t - b):
                annot.xytext = (annot.xytext[0], -annot.xytext[1])
            """
            text = "ID:  {}\nSEQ: {}\nΔv:  {}".format(idList[idx], seq[idx], round(Δv[idx],2))
            annot.set_text(text)
            # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
            annot.get_bbox_patch().set_alpha(0.4)


        def hover(event):
            vis = annot.get_visible()
            if event.inaxes == ax:
                cont, ind = sc.contains(event)
                if cont:
                    update_annot(ind)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

        fig.canvas.mpl_connect("motion_notify_event", hover)

        plt.show()

    def export(self, filename = 'auto', exportFullTree = False):
        # IMPORTINT PACKAGES
        pd = self.pd
        spk = self.spk

        # FINDING RESUTLS
        id_ = self.getResults(printR = False)

        rowCID = []
        for id in id_:
            lineage = self.getLineage(id)[:-1]
            for i in lineage:
                if i not in rowCID:
                    rowCID.append(i)

        # print(rowCID)

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
        out = {'id': []}
        with self.click.progressbar(id_, fill_char = '█', empty_char = ' ', label='Exporting Results    ') as bar:
            for id in bar:
                out['id'].append(id)
                attr = self.queryBranchAttributes(id)

                for key in attr.keys():
                    if key in out.keys():
                        out[key].append(attr[key])
                    else:
                        out[key] = []
                        out[key].append(attr[key])

        # CONVERTING TO PANDAS & SAVING
        try:
            df = pd.DataFrame(out)
        except Exception as inst:
            print(inst)
            return out

        df.to_excel(writer, sheet_name = 'Results')

        # EXPORTING FULL TREE
        if exportFullTree:
            thresh = 1000000
            out = {'id': [], 'layer': [], 'visits': [], 'terminal': [], 'planet': [], 'date': [], 'theta': [], 'delta': [], 'h': [], 'dvTot': [], 'dvLeg': [], 'x': [], 'parent': [], 'children': [], 'vinf': []}
            with self.click.progressbar(range(len(self.node)), fill_char = '█', empty_char = ' ', label='Exporting All Nodes  ') as bar:
                for i in bar:
                    if i == thresh-1:
                        # EXPORTING TO SPREADSHEET
                        df = pd.DataFrame(out)

                        sheetName = "All Nodes (1)"
                        df.to_excel(writer, sheet_name = sheetName)

                        for key in out.keys():
                            out[key] = []

                    else:
                        # NODE INFORMATION
                        out['id'].append(i)
                        out['layer'].append(self.node[i].layer)
                        out['visits'].append(self.node[i].n)
                        out['terminal'].append("True" if self.node[i].isTerminal else "False")

                        # NODE STATE
                        if self.node[i].state and len(self.node[i].state) == 2:
                            out['planet'].append(self.pkP[self.node[i].state[0]])
                            out['date'].append(spk.et2utc(self.node[i].state[1], 'C', 14, 12))
                            out['theta'].append("None")
                        elif self.node[i].state:
                            K = self.node[i].state[0][1]
                            p = self.pkP[self.node[i].state[0][0]]
                            out['planet'].append("{} {} DSM".format(p, K))
                            out['date'].append(spk.et2utc(self.node[i].state[1], 'C', 14, 12))
                            out['theta'].append(self.node[i].state[2])
                        else:
                            out['planet'].append("None")
                            out['date'].append("None")
                            out['theta'].append("None")

                        # NODE COST/REWARDS
                        out['dvTot'].append(round(self.node[i].Δv, 2))
                        out['dvLeg'].append(round(self.node[i].ΔvLeg, 2))
                        out['delta'].append(round(self.node[i].delta, 2) if self.node[i].delta else [])
                        out['h'].append(round(self.node[i].h.item(), 2) if self.node[i].h and self.node[i].h != -1 else [])
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
                            vp = self.spk.spkezr(self.node[i].state[0][0], self.node[i].state[1], self.frame, self.abcorr, self.cntrBody)[0][3:]
                            vinf = self.np.linalg.norm(self.np.array(l.get_v2()[0])/1000 - vp)
                            out['vinf'].append(vinf)
                        else:
                            out['vinf'].append([])


                # EXPORTING TO SPREADSHEET
                df = pd.DataFrame(out)
                """
                def colorRow(x):
                    print("T")
                    df2 = x.copy()
                    for i in rowCID:
                        df2.loc[df2['id'] == rowCID[i], :] = 'background-color: yellow'
                    return df2

                # try:
                print("APLYING STYLE")
                df.style.apply(colorRow, axis = None)
                # except:
                    # print("jinja2 Package Not Found -- No styling will be applied")
                """
                sheetName = "All Nodes" if len(self.node) < thresh-1 else "All Nodes (2)"
                df.to_excel(writer, sheet_name = sheetName)

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

    def exit(self, *args):
        import signal
        import sys
        sysExit = lambda *args: sys.exit(1)
        signal.signal(signal.SIGTSTP, sysExit)
        try:
            r = input("\n\nAre you sure you want to exit, all unsaved information will be lost? [Yes/No] ")
            if r == "Yes":
                sys.exit(1)
        except KeyboardInterrupt:
            print("\nForce Quitting")
            sys.exit(1)

    def debug(self, *args):
        print(" ")
        print("### =================================================== ###")
        print("### ==================== DEBUGGING ==================== ###")
        print("### =================================================== ###")
        print(" ")
        self.breakLoop = True
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
        self.delta = []
        # Establishes Parent above (Single Scalar of IDs)
        self.parent = parent
        self.children = None  # Children Below (Array of IDs)
        self.state = state  # Either NAIF ID String or Epoch
        self.layer = layer
        self.isTerminal = isTerminal