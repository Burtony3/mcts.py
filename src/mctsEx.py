# IMPORTING TOOLBOX
from MCTS import MCTS

# INITIALIZING OBJECT
mcts = MCTS()

# SETTING UP CONSTRAINTS
mcts.setup(
    '7',
    ["Jan 01, 2037", "Jan 02, 2038"],
    maxIters = 1000,
    dvBudget = 10,
    detail = 8,
    flybyBodies = ['3', '5'],
    C3max = 45,
    minimizeVinfCap = True,
    freeCapture = False
)

# RUNNING SIMULATION
mcts.run()

# PRINTING RESULTS
idList = mcts.getResults()

# SAVING RESULTS
mcts.export(filename = 'out.xlsx', exportFullTree = True)

# PLOTTING BEST
mcts.plotFamilies(idList, hideAxes = False, showTop = 15)

mcts.plotPath(idList[0])
