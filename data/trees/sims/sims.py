# IMPORTING TOOLBOX
from MCTS import MCTS

# INITIALIZING OBJECT
mcts = MCTS()

# SETTING UP CONSTRAINTS
mcts.setup(
    '6',
    ["Jan 01, 2020", "Jan 01, 2023"],
    maxIters = 5000,
    dvBudget = 10,
    detail = 16,
    C3max = 45,
    enableDSM = True,
    flybyBodies = ['2'],
)

# RUNNING SIMULATION
mcts.run()

# PRINTING RESULTS
idList = mcts.getResults()

# print(mcts.queryBranchAttributes(idList[0]))

# SAVING RESULTS
mcts.export(filename = 'sims.xlsx', exportFullTree = True)

# PLOTTING BEST
mcts.plotFamilies(idList, hideAxes = False, showTop = 20)

# mcts.plotPath(idList[0])
mcts.plotAttributes(idList, 'dateL', 'seq', 'C3')
