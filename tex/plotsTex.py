import spiceypy as spk
import pykep as pk
from pykep import propagate_lagrangian
import matplotlib
from matplotlib import patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import dill
import sys
sys.path.insert(0, "/home/burtonyale/Documents/repos/MCTS.py/src")
f = open("../data/trees/clipper/clipper70k.pckl", "rb")
tree = dill.load(f)
tree.loadKernels()

def state(tree, xAttr, yAttr, cAttr, showTop = 50, cmap = "turbo", action = "show", removeVinfdv = True):
    # plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    # plt.rc('text', usetex=True)
    # ADJUSTING LABELS FOR LATEX
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif', size = 14)
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{xfrac}"

    # SETTING DICTIONARIES
    # (Label, Grid Bool)
    lDict = {
        "C3":    (r"Launch C3 ($\displaystyle\sfrac{km^2}{s^2}$)", "on"),
        "dateL": (r'Launch Date', "on"),
        "tof":   (r'Total Time of Flight (Days)', "on"),
        "dvTot": (r"Unoptimized Mission $\displaystyle\Delta$v ($\displaystyle\sfrac{km}{s}$)", "on"),
        "vInfF": (r"Arrival $\displaystyle V_\infty$ ($\displaystyle\sfrac{km}{s}$)", "on"),
        "seq":   (r'Planetary Arrival Sequence', "off"),
        }
    shapeList = ["o", "^", "s", "*", "X", "d", "h", ".", "1", "+", "p"]
    cList = ["blue", "orange", "green", "red", "purple", "gray", "pink", "cyan", "olive", "brown"]
    cList = ["tab:"+color for color in cList]

    # STARTING PLOT WINDOW
    fig = plt.figure()
    ax = fig.gca()

    # SETTING LABELS & TITLES
    ax.set_xlabel(lDict[xAttr][0])
    ax.set_ylabel(lDict[yAttr][0])


    # INITIALIZING PLOTTED VALUES
    x = []
    y = []
    c = []
    seq = []
    Δv = []

    # RETRIEVING PLOTTED VALUES
    idList = tree.getResults()
    if type(showTop) == int:
        idList = idList[showTop::-1]
    else:
        idList = idList[::-1]
    for id in idList:
        attr = tree.queryBranchAttributes(id)
        x.append(attr[xAttr])
        if yAttr == 'dvTot' and removeVinfdv:
            y.append(attr[yAttr] - attr['vInfF'])
        else:
            y.append(attr[yAttr])
        if cAttr and cAttr[:3] != 'top':
            c.append(attr[cAttr])
        seq.append(attr['seq'])
        if removeVinfdv:
            Δv.append(attr['dvTot'] - attr['vInfF'])
        else:
            Δv.append(attr['dvTot'])

    # PLOTTING
    ax.grid(lDict[xAttr][1], 'major', 'x', linestyle = "--", zorder = 3)
    ax.grid(lDict[yAttr][1], 'major', 'y', linestyle = "--", zorder = 3)
    sDict = {}
    cDict = {}
    itr = 0
    if cAttr == "seq" or cAttr[:3] == 'top': 
        for i in range(len(x)):
            newBool = False
            # CHECKING FOR COLOR MATCHES KEY
            if cAttr == 'seq':
                size = 20
                if seq[i] in sDict.keys():
                    shape = sDict[seq[i]]
                    color = cDict[seq[i]]
                else:
                    shape = shapeList[itr]
                    color = cList[itr]
                    sDict[seq[i]] = shape
                    cDict[seq[i]] = color
                    itr += 1
                    newBool = True
            else:
                size = 5
                shape = 'o'
                color = 'gray' if i > len(x) - int(cAttr[3:]) else 'black'
            if newBool:
                sc = plt.scatter(x[i], y[i], s=size, zorder = 3, marker = shape, label = seq[i], color = color)
            else:
                sc = plt.scatter(x[i], y[i], s=size, zorder = 3, marker = shape, color = color)
    else:
        if cAttr == None:
            size = 5
            sc = plt.scatter(x, y, s=size, zorder = 3, color = "black")
        else:
            size = 7.5
            sc = plt.scatter(x, y, c = c, s=size, zorder = 3, cmap=cmap)
            cbar = fig.colorbar(sc)
            cbar.ax.set_ylabel(lDict[cAttr][0])
    # print(sDict)

    if xAttr == "seq" or xAttr == "dateL":
        # print(vals, labels)
        if xAttr == "dateL": 
            vals, labels = plt.xticks()
            labels = [spk.et2utc(val, 'C', 14, 12) for val in vals]
            plt.xticks(vals, labels, rotation = 30, ha = 'right', size=12)
        plt.xticks(rotation = 30, ha = 'right', size=12)

    if cAttr == "seq":
        plt.legend(
            numpoints = 1, 
            loc = "lower left", 
            ncol = len(sDict) // 2 + (1 if len(sDict) % 2 == 1 else 0), 
            handletextpad = 0.1, prop={'size': 12}, 
            bbox_to_anchor = (0, 1.02, 1, 0.2)
        )

    # ADDING COLORBAR
    plt.tight_layout()

    if action == "show":
        plt.show()
    else:

        plt.savefig("filename.png")


def orbit(tree, showTop = 50, N = 60, id = None, seqT = None):
    AU2m = 149597870700
    AU2km = AU2m/1000
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif', size = 14)
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{xfrac}"
    xLabel = "X-Position J2000 Frame (AU)"
    yLabel = "Y-Position J2000 Frame (AU)"

    # STARTING PLOT WINDOW
    fig = plt.figure()
    ax = fig.gca()

    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)

    # CREATING COLOR LISTS
    cList = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan"]
    cList = ["tab:"+color for color in cList]
    cDict = {}

    if not id:
        idList = tree.getResults()
        if type(showTop) == int:
            idList = idList[showTop::-1]
        else:
            idList = idList[::-1]
        id = idList[0]
    else:
        idList = id

    # FINDING INITIAL EPOCH
    et0 = tree.node[id].state[1]
    if tree.enableDSM: tree.P.append("3") # Adding Earth for DSM Case
    for P in tree.P:
        # FINDING SECONDS IN 1 ORBIT PERIOD
        tau = tree.tau[P]*86400
        et1 = et0 + tau
        et = tree.np.linspace(et0, et1, 250)

        # COLLECTING DATA POINTS
        x = []
        y = []
        for et_ in et:
            r = tree.spk.spkezr(P, et_, tree.frame, tree.abcorr, tree.cntrBody)[0][:3]/tree.AU2km
            x.append(r[0])
            y.append(r[1])
        ax.plot(x, y, ':k', linewidth=1)

    # PLOTTING CENTER BODY SURFACE
    ax.scatter(0, 0, s = 50, c = "#FFCC33", zorder = 4)

    # RETRIEVING PLOTTED VALUES
    itr = 0
    for id in idList:
        attr = tree.queryBranchAttributes(id)
        seq = attr['seq']
        lineage = tree.getLineage(id)[-3::-1]

        if not seqT or seq in seqT:
            if seq in cDict.keys():
                color = cDict[seq]
            else:
                color = cList[itr]
                cDict[seq] = color
                itr += 1

            for id_ in lineage:
                l = tree.computeLambert(id_)
                r = l.get_r1()
                v = l.get_v1()[0]
                T = l.get_tof()
                mu = l.get_mu()

                dt = T / (N - 1)

                x = np.array([0.0] * N)
                y = np.array([0.0] * N)

                for i in range(N):
                    x[i] = r[0] / AU2m
                    y[i] = r[1] / AU2m
                    # print(r, type(r), v, type(v), dt, type(dt), mu, type(mu))
                    r, v = propagate_lagrangian(r, v, dt, mu)

                ax.plot(x, y, linewidth = 1, color = color)

    # print(cDict)

    # ADDING FAMILIES TO LEGENDS
    patchList = []
    for key in cDict:
        # lab = []
        # for i in range(len(key)):
        #     if key[i] in tmpDict.values():
        #         lab.append([key_ for key_ in tmpDict.items() if key_[1] == key[i]][0][0])
        #     else:
        #         lab.append(self.letterDict[key[i]] if key[i] in self.letterDict else key[i])
        # lab = ''.join(lab)
        dataKey = mpatches.Patch(color = cDict[key], label=key)
        patchList.append(dataKey)

    if len(patchList) > 1:
        plt.legend(handles = patchList)

    ax.grid(zorder = 3)
    ax.set_aspect('equal', 'box')
    # ADDING COLORBAR
    plt.tight_layout()
    plt.show()


# state(tree, 'seq', 'dvTot', "top50", showTop = None)
orbit(tree, showTop = 25, seqT = ["EEVEEJ"])