class MCTS(object):
    def __init__(maxIter, dvBudget, maxNumFlyby, launchWindow, finalBody; kwargs...):
        # Input Parse


        # Creating Root Node
        self.node = []              # Allocates List
        self.node[0] = nodeNew('3') # Root of all Nodes (nodeNew will replace node.py)

        # Running to Completion
        for i in range(maxIters):
            # Select
            id, expandBool = select()   # Loops until finds leaf node & boolean if expand

            # Expand?
            if expandBool is True:
                id = expand(id, finalBody)

            # Simulate
            value = simulate(id, finalBody, maxNumFlyby, dvBudget)

            # Backprop
            backprop(id, value)

        # Returning Results
        return node
