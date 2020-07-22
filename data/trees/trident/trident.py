# IMPORTING TOOLBOX
from MCTS import MCTS

# INITIALIZING OBJECT
mcts = MCTS()

# SETTING UP CONSTRAINTS
mcts.setup(
    '8',
    ["Jan 01, 2029", "Jul 01, 2029"],
    maxIters = 10000,
    dvBudget = 4.5,
    detail = 16,
    flybyBodies = ['2', '5'],
    C3max = 35,
    enableDSM = True,
    freeCapture = True,
)

# RUNNING SIMULATION
mcts.run()

# PRINTING RESULTS
idList = mcts.getResults()

# print(mcts.queryBranchAttributes(idList[0]))

# SAVING RESULTS
mcts.export(filename = 'trident.xlsx', exportFullTree = True)

# PLOTTING BEST
mcts.plotFamilies(idList, hideAxes = False, showTop = 20)

# mcts.plotPath(idList[0])
mcts.plotAttributes(idList, 'dateL', 'tof', 'dvTot')
