def getLineage(id):
    lineage = []
    while id is not 0:
        lineage.append(id)
        id = node[id].parent
    lineage.append(0)

    # Returns most recent node @ idx = 0 & oldest parent @ idx = -1
    return lineage
