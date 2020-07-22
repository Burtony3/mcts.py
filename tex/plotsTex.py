import spiceypy as spk
import pykep as pk
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import dill
import sys
sys.path.insert(0, "/home/burtonyale/Documents/repos/MCTS.py/src")
f = open("../data/trees/clipper/clipper60k.pckl", "rb")
tree = dill.load(f)
tree.loadKernels()

def state(tree, xAttr, yAttr, cAttr, showTop = 50, cmap = "turbo", action = "show"):
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
        y.append(attr[yAttr])
        if cAttr:
            c.append(attr[cAttr])
        seq.append(attr['seq'])
        Δv.append(attr['dvTot'])

    # PLOTTING
    ax.grid(lDict[xAttr][1], 'major', 'x', linestyle = "--", zorder = 3)
    ax.grid(lDict[yAttr][1], 'major', 'y', linestyle = "--", zorder = 3)
    sDict = {}
    cDict = {}
    itr = 0
    if cAttr == "seq": 
        size = 20
        for i in range(len(x)):
            newBool = False
            # CHECKING FOR COLOR MATCHES KEY
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


state(tree, 'dateL', 'dvTot', "seq", showTop = None)