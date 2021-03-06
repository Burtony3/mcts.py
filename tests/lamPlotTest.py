# MATH TOOLBOXES
import numpy as np
import spiceypy as spk   # cSpice Kernels
from scipy.optimize import newton
import pykep as pk
import math

# PLOTTING TOOLBOXES
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# ===================== USER INPUTS ===================== #
Δ    = 24               # Angular Resolution (Δθ = 360/Δ)
p1   = '3'              # Origin Planet NAIF ID
utc1 = "Oct 21, 2024"   # Launch Date
p2   = '3'              # Flyby Planet NAIF ID
p3   = '5'              # Arrival Planet NAIF ID
# ======================================================= #

# SPICE SETUP
spk.furnsh("../data/spk/de438.bsp")
spk.furnsh("../data/spk/naif0009.tls")
spk.furnsh("../data/spk/gm_de431.tpc")
spk.furnsh("../data/spk/pck00010.tpc")
frame = "ECLIPJ2000"
abcorr = "NONE"
J2000_jd = 2451544.5

# DICTIONARIES
pkP = {'1': 'mercury', '2': 'venus',   '3': 'earth', '4': 'mars', '5': 'jupiter', '6': 'saturn',  '7': 'uranus',  '8': 'neptune'}
tau = {'1': 87.97,     '2': 224.7,     '3': 365.25,  '4': 687,    '5': 4331,      '6': 10747,     '7': 30589,     '8': 59800}
col = {'1': '#5a6165', '2': '#955c1c', '3': 'b',     '4': 'r',    '5': '#9e7a5b', '6': '#c9c0a9', '7': '#8eb0b8', '8': '#4d80d7'}

# 3D PLOT AXIS EQUIL
def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

# UPDATING FIGURE ELEMENTS
def updateLam():
    # IMPORTING SLIDERS AND ADJUSTING INPUTS
    global idx1
    global idx2
    ax.cla()
    l_ = l[idx1]
    l2_ = l2[idx1][idx2]
    et2_ = et2[idx1]
    et3_ = et3[idx1][idx2]

    # PLOTTING ORIGIN BODY
    pl = pk.planet.jpl_lp(pkP[p1])
    pk.orbit_plots.plot_planet(
        pl, 
        t0 = pk.epoch(float(spk.et2utc(et1, 'J', 10)[3:]) - J2000_jd), 
        axes = ax, 
        color=col[p1]
    )

    # PLOTTING FLYBY BODY
    if int(p2) < 5:
        pl = pk.planet.jpl_lp(pkP[p2])
        pk.orbit_plots.plot_planet(
            pl, 
            t0 = pk.epoch(float(spk.et2utc(et2_, 'J', 10)[3:]) - J2000_jd), 
            axes = ax, 
            color=col[p2]
        )

    # PLOTTING ARRIVAL BODY
    if int(p3) < 5:
        pl = pk.planet.jpl_lp(pkP[p3])
        pk.orbit_plots.plot_planet(
            pl, 
            t0 = pk.epoch(float(spk.et2utc(et3_, 'J', 10)[3:]) - J2000_jd), 
            axes = ax, 
            color=col[p3]
        )

    # PLOTTING LAMBERT ARC
    pk.orbit_plots.plot_lambert(l_, axes = ax, color='c')
    pk.orbit_plots.plot_lambert(l2_, axes = ax, color='c')

    # GETTING TRAJECTORY PROPERTIES
    C3, tof, Δv = getLamProps(l_, l2_, s1, s2[idx1])

    # APPLYING TO FIGURE
    # ax.set_title("LAUNCH C3: {:.2f} $km^2/s^2$\nTIME OF FLIGHT: {:.2f} + {:.2f} DAYS\nCORRECTION Δv: {:.2f} $km/s$".format(C3, tof[0], tof[1], Δv))
    d2 = spk.et2utc(et2_, 'C', 14, 12)
    d3 = spk.et2utc(et3_, 'C', 14, 12)
    strC3 = "LAUNCH C3: {:.2f} $km^2/s^2$\n".format(C3)
    strFlyby = "FLYBY DATE: {:s}\n".format(d2)
    strArr = "ARRIVAL DATE: {:s}\n".format(d3)
    strDV = "REQ. FLYBY Δv: {:.2f} $km/s$".format(Δv)
    ax.set_title(strC3 + strFlyby + strArr + strDV)
    fig.canvas.draw_idle()
    axisEqual3D(ax)

# GETS LAMBERT ARC PROPERTIES
def getLamProps(l, l2, state1, state2):
    """
    NOTE: All pykep lambert calculations are done in meters & seconds and are converted to km
    """
    C3 = np.linalg.norm(np.array(l.get_v1()[0])/1000 - state1[3:6])**2 # Launch C3 (km²/s²)
    tof = [l.get_tof()/86400, l2.get_tof()/86400]                      # Time of Flight (days)
    vih = np.array(l.get_v2()[0])/1000
    voh = np.array(l2.get_v1()[0])/1000
    vi = vih - state2[3:6]                    # Planet 2 Arc 1 v∞
    vo = voh - state2[3:6]                   # Planet 2 Arc 2 v∞

    mu = spk.bodvrd(p2 + "99", "GM", 1)[1]
    aOutI = -mu / np.linalg.norm(vo)**2

    def f(eOut):
        aOut = -mu/np.linalg.norm(vo)**2
        aIn = -mu/np.linalg.norm(vi)**2
        delta =  math.acos(np.dot(vo, vi) / (np.linalg.norm(vo)*np.linalg.norm(vi)))
        # print(aOut, aIn, eOut, delta)
        eOut = (aOut / aIn) * (eOut - 1) * math.sin(delta - math.asin(1 / eOut)) - 1
        return eOut


    def f_prime(eOut):
        aOut = -mu/np.linalg.norm(vo)**2
        aIn = -mu/np.linalg.norm(vi)**2
        delta = math.acos(np.dot(voh, vi) / (np.linalg.norm(voh)*np.linalg.norm(vi)))
        tmp = [aOut / aIn, delta - math.asin(1 / eOut)]
        eOut = (tmp[0] * (eOut - 1) + 1) * (math.cos(tmp[1]) / (eOut**2 * math.sqrt(1 - eOut**-2))) + tmp[0] * math.sin(tmp[1])
        return eOut

    try:
        rp = aOutI * (1 - newton(f, 1.5, f_prime))
        print(rp - spk.bodvrd(p2 + "99", "RADII", 3)[1][1])
        tmp = (2 * mu) / rp
        norm = np.linalg.norm
        Δv = norm(math.sqrt(norm(vo)**2 + tmp)) - norm(math.sqrt(norm(vi)**2 + tmp))
    except:
        Δv = 999.0
    # Δv = np.linalg.norm(vo - vi)                                       # Norm difference in above 2
    return C3, tof, abs(Δv)

def setEpoch(t0, p0, p1):

    if p0 == p1:      # If returning to same planet
        det = Δ // 3
        mod = 3 - Δ % 3
        et = np.array([])
        i = 1
        for K in range(1, 4):
            det_ = det + 1 if i > mod else det
            i += 1
            et0 = 0.9  * K * tau[p0[0]]*86400 + t0
            et1 = (K * tau[p0[0]] - 3)*86400 + t0
            et_ = np.linspace(et0, et1, det_)
            for val in et_: et = np.append(et, val)
    
    else:
        if int(p1) > 4:   # Outer Planets Case
            n = 0.05
            m = 0.33
        else:             # All other cases
            n = 0.1
            m = 1

        # SETTING EPOCH UPPER AND LOWER LIMITS
        et0 = n*(tau[p0[0]] + tau[p1[0]])*86400 + t0
        et1 = m*(tau[p0[0]] + tau[p1[0]] - 1)*86400 + t0

        # CREATING LINEAR RANGE
        et = np.linspace(et0, et1, Δ)
        """
    # SETTING PARAMETERS FOR RANGES
    if p0 == p1: # If returning to same planet
        n = 0.495 # Return in ~2.1 periods
        m = 1.5  # Return in 3 periods
    else:        # All other cases
        n = 0.1
        m = 1

    # SETTING EPOCH UPPER AND LOWER LIMITS
    et0 = n*(tau[p0] + tau[p1])*86400 + t0
    et1 = m*(tau[p0] + tau[p1] - 1)*86400 + t0

    # CREATING LINEAR RANGE
    et = np.linspace(et0, et1, Δ)
        """
    return et

# EVALUATES ARC #1 SLIDER CHANGE
def updateA(val):
    global idx1
    idx1 = val
    updateLam()

# EVALUATES ARC #2 SLIDER CHANGE
def updateB(val):
    global idx2
    idx2 = val
    updateLam()

# ORIGIN BODY STATE
et1 = spk.str2et(utc1)
s1 = spk.spkezr(p1, et1, frame, abcorr, '0')[0]

# FLYBY BODY STATES (size = [Δ, 1])
et2 = setEpoch(et1, p1, p2)
s2 = [spk.spkezr(p2, et2[i], frame, abcorr, '0')[0] for i in range(Δ)]

# ARRIVAL BODY STATES (size = [Δ, Δ])
et3 = [] # Starting List
s3 = []  # Starting List
for i in range(Δ):      # Looping through all flyby states
    et3.append(setEpoch(et2[i], p2, p3))
    for j in range(Δ):  # Looping through all arrival states
        s3.append(tuple(spk.spkezr(p3, et3[i][j], frame, abcorr, '0')[0]))
s3 = [s3[Δ*i:Δ*i+Δ] for i in range(Δ)]   # Converting to 2D List

# COMPUTING LAMBERTS
l = []
l2 = []
for i in range(Δ):
    l.append(pk.lambert_problem(
        r1 = s1[0:3]*1000,
        r2 = s2[i][0:3]*1000,
        tof = et2[i]-et1,
        mu = pk.MU_SUN
    ))
    for j in range(Δ):
        tof_ =  et3[i][j]-et2[i]
        l2.append(pk.lambert_problem(
            r1 = s2[i][0:3]*1000,
            r2 = np.asarray(s3[i][j][0:3])*1000,
            tof = tof_,
            mu = pk.MU_SUN
        ))

l2 = [l2[Δ*i:Δ*i+Δ] for i in range(Δ)]



idx1 = 0
idx2 = 0
# idx1 = 37
# idx2 = 17

# ASSIGNING INITIAL LAMBERT PLOT INDEX
l_ = l[idx1]
l2_ = l2[idx1][idx2]
et2_ = et2[idx1]
et3_ = et3[idx1][idx2]

# STARTING PLOT WINDOW
fig = plt.figure()
ax = fig.gca(projection = '3d')
ax.view_init(azim=0, elev=90)

# PLOTTING FIRST LAMBERT
updateLam()

# CREATING SLIDER OBJECTS
ax_a = plt.axes([0.25, 0.08, 0.5, 0.025])
ax_b = plt.axes([0.25, 0.05, 0.5, 0.025])
sl_1 = Slider(ax_b, 'Arc 1 #', 0, len(l)-1, valstep=1, valinit=idx1, valfmt='%i')
sl_2 = Slider(ax_a, 'Arc 2 #', 0, len(l)-1, valstep=1, valinit=idx2, valfmt='%i')

# ASSINGING VALUE CHANGED CALLBACK
sl_1.on_changed(updateA)
sl_2.on_changed(updateB)

# SHOWING PLOT
plt.show()