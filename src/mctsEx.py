# IMPORTING TOOLBOX
from MCTS import MCTS

# INITIALIZING OBJECT
mcts = MCTS()

# SETTING UP CONSTRAINTS
mcts.setup(
    '5',
    ["Jan 01, 2037", "Dec 31, 2037"],
    maxIters = 10000,
    dvBudget = 10,
    detail = 16,
    flybyBodies = ['2', '3'],
    C3max = 30,
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
mcts.plotFamilies(idList, hideAxes = False, showTop = 50)
