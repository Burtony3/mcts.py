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
f = open("../data/trees/trident/triton.pckl", "rb")
tree = dill.load(f)
tree.loadKernels()
# idList = np.loadtxt("idList.txt", delimiter="\n")


def state(tree, xAttr, yAttr, cAttr, showTop = 50, cmap = "turbo", action = "show", removeVinfdv = True, pickID = None):
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
        if pickID and id == pickID:
            pIDx = x[-1]
            pIDy = y[-1]

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

    if pickID:
        ax.scatter(pIDx, pIDy, s = 65, edgecolors = "black", marker = 'o', facecolors = 'none', linewidths = 1.25, label = "Plotted Result", zorder = 5)

    if xAttr == "seq" or xAttr == "dateL":
        # print(vals, labels)
        if xAttr == "dateL": 
            vals, labels = plt.xticks()
            labels = [spk.et2utc(val, 'C', 14, 12) for val in vals]
            plt.xticks(vals, labels, rotation = 30, ha = 'right', size=12)
        plt.xticks(rotation = 30, ha = 'right', size=12)

    if cAttr == "seq":
        ncol = len(sDict) // 2 + (1 if len(sDict) % 2 == 1 else 0)
        if len(sDict) % 2 == 0 and pickID != None: ncol += 1
        plt.legend(
            numpoints = 1, 
            loc = "lower left", 
            ncol = ncol, 
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
    def plotLam(l):
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
        return x, y

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

    if type(id) != np.ndarray and type(id) != list:
        idList = tree.getResults()
        if type(showTop) == int:
            idList = idList[showTop::-1]
        else:
            idList = idList[::-1]
        id = idList[0]
    else:
        idList = [int(val) for val in id]
        id = idList[0]
        if len(idList) > showTop:
            idList = idList[showTop::-1]
        else:
            idList = idList[::-1]


    # FINDING INITIAL EPOCH
    print(idList)
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
                if len(tree.node[id_].state) == 3:
                    # RETRIEVING PARENT AND STATE
                    parent = tree.node[id_].parent
                    X = tree.spk.spkezr(tree.node[parent].state[0][0], tree.node[parent].state[1], tree.frame, tree.abcorr, tree.cntrBody)[0]

                    # RETRIVING DSM INFORMATION
                    K = tree.node[id_].state[0][1]
                    _, _, dt, vinfL = tree.queryDSM(K, tree.node[id_].state[2])

                    # CALCULATING DSM APEHELION
                    rp1 = X[:3]                         # Perihelion position
                    theta = tree.getJ2000Ang(rp1)       # Angle of perihelion wrt J2000
                    vp = X[3:]                          # Perihleion velocity
                    r = tree.np.linalg.norm(rp1)        # Normal Value
                    v = tree.np.linalg.norm(vp) + vinfL # Normal Value
                    mu = tree.spk.bodvrd('sun', 'GM', 1)[1].item() # μ of Sun
                    a = -mu/(v**2 - (2*mu)/r)           # Semi-Major Axis
                    tof = tree.math.pi*tree.math.sqrt(a**3 / mu) # ½ period time of flight
                    ra = 2*a - r                        # Apehelion normal value           
                    phi = theta + tree.math.pi + tree.np.radians(1) # Angle of apehelion wrt J2000 (+ small bias)
                    ra = [ra*tree.math.cos(phi), ra*tree.math.sin(phi), 0]  # Apehelion vector

                    # MARKING DSM MANEUVER POINT
                    ax.scatter(ra[0]/tree.AU2km, ra[1]/tree.AU2km, color = color, marker='x')

                    # CREATING LAMBERT ARC INPUTS
                    ra = tree.np.array(ra)
                    rp2 = tree.spk.spkezr(tree.node[id_].state[0][0], tree.node[id_].state[1], tree.frame, tree.abcorr, tree.cntrBody)[0][:3]
                    
                    # LAMBERT ARC FROM LAUNCH TO DSM POINT
                    l1 = tree.pk.lambert_problem(
                        r1 = rp1*1000,
                        r2 = ra*1000,
                        tof = tof,
                        mu = tree.pk.MU_SUN
                    )
                    x, y = plotLam(l1)
                    ax.plot(x, y, linewidth = 1, color = color)

                    # LAMBERT ARC FROM DSM POINT TO FLYBY
                    l2 = tree.pk.lambert_problem(
                        r1 = ra*1000,
                        r2 = rp2*1000,
                        tof = (dt*86400 - tof),
                        mu = tree.pk.MU_SUN
                    )
                    x, y = plotLam(l2)
                    ax.plot(x, y, linewidth = 1, color = color)

                else:
                    l = tree.computeLambert(id_)

                    x, y = plotLam(l)

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
        #         lab.append(tree.letterDict[key[i]] if key[i] in tree.letterDict else key[i])
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


# state(tree, 'tof', 'dvTot', "seq", showTop = 100, pickID = 384842)
orbit(tree, id = [121249, 109667])