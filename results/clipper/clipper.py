# IMPORTING TOOLBOX
from MCTS import MCTS

# INITIALIZING OBJECT
mcts = MCTS()

# SETTING UP CONSTRAINTS
mcts.setup(
    '5',
    ["Mar 01, 2022", "Sept 01, 2022"],
    maxIters = 75000,
    dvBudget = 15,
    detail = 24,
    C3max = 10,
    enableDSM = False,
)

# RUNNING SIMULATION
mcts.run()

# PRINTING RESULTS
idList = mcts.getResults()

# print(mcts.queryBranchAttributes(idList[0]))

# SAVING RESULTS
mcts.export(filename = "clipper.xlsx", exportFullTree = True)

# PLOTTING BEST
mcts.plotFamilies(idList, hideAxes = False, showTop = 20)
