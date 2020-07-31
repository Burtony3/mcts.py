# IMPORTING TOOLBOX
from MCTS import MCTS

# INITIALIZING OBJECT
mcts = MCTS()

# SETTING UP CONSTRAINTS
mcts.setup(
    '5',
    ["Jun 01, 1989", "Dec 31, 1991"],
    maxIters = 50000,
    dvBudget = 10,
    detail = 16,
    C3max = 20,
    enableDSM = True,
)

# RUNNING SIMULATION
mcts.run()

# PRINTING RESULTS
idList = mcts.getResults()

# print(mcts.queryBranchAttributes(idList[0]))

# SAVING RESULTS
mcts.export(filename = 'galileo.xlsx', exportFullTree = True)

# PLOTTING BEST
mcts.plotFamilies(idList, hideAxes = False, showTop = 20)

# mcts.plotPath(idList[0])
mcts.plotAttributes(idList, 'dateL', 'seq', 'C3')
