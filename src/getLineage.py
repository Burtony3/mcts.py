def getLineage(id):
    lineage = []
    while id is not 0:
        lineage.append(id)
        id = node[id].parent
    lineage.append(0)
    return lineage
