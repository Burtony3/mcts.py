import math

# Selection
def select():
    # Constant Initialization
    cp = 1/math.sqrt(2.0) # Cost Adjustment Parameter
    id = 0                # Starting at top of tree
    node = self.node      # FIXME: Redefining node (Necessary?)

    while node[id].children != None:
        ucb1 = []
        
        # Checks Node's children
        for i in range(len(node[id].children)):
            X = node[node[id].children[i]].cost
            N = node[id].n
            n = node[node[id].children[i]].n

            if n == 0: # Immediatly selects first unvisited node
                id = node[id].children[i]
                break
            else:
                ucb1.append(X + cp*math.sqrt(math.log1p(N)/n))

        # Find UCB w/ max value
        indexChild = ucb1.index(max(ucb1))

        # Set the id to this child's index in the Node list
        id = node[id].children[indexChild]

    if node[id].n == 0:
        expandBool = False
    else:
        expandBool = True

    return id , expandBool
