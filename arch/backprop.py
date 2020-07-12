def backprop(id, value):

    lineage = getLineage(id)

    for id_ in lineage:
        node[id_].cost += value
        node[id_].n += 1
