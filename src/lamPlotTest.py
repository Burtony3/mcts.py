import sys
import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spk   # cSpice Kernels
import pykep as pk
from matplotlib.widgets import Slider
        

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

# SPICE SETUP
spk.furnsh("../data/spk/de438.bsp")
spk.furnsh("../data/spk/naif0009.tls")
frame = "ECLIPJ2000"
abcorr = "NONE"
J2000_mjd = 51544

# EARTH DEPARTURE
etE = spk.str2et("Jun 01, 2020")
mdj2000E = float(spk.et2utc(etE, 'J', 10)[3:]) - 2451544.5
# print(pk.epoch(mdj2000E))
stateE = spk.spkezr('3', etE, frame, abcorr, '0')[0]

# MARS FLYBY
tauM = 687
Δ = 16
etM = [etE + (i+1)*(687*86400/Δ) for i in range(Δ)]
stateM = [spk.spkezr('4', etM[i], frame, abcorr, '0')[0] for i in range(Δ)]

# JUPITER ARRIVAL
etJ = np.zeros([Δ, Δ])
stateJ = []
for i in range(len(etM)):
    for j in range(Δ):
        etJ[i, j] = etM[i] + (j+1)*(4331*86400/Δ)
        stateJ.append(tuple(spk.spkezr('5', etJ[i, j], frame, abcorr, '0')[0]))
# print(len(stateJ), type(stateJ), type(stateJ[0]))
# print(etJ[:32])
# print(etJ[:2])
stateJ = [stateJ[i:i+Δ] for i in range(Δ)]

# COMPUTING LAMBERTS
l = []
l2 = []
for i in range(Δ):
    l.append(pk.lambert_problem(
        r1 = stateE[0:3]*1000,
        r2 = stateM[i][0:3]*1000,
        tof = etM[i]-etE,
        mu = pk.MU_SUN
    ))
    for j in range(Δ):
        tof_ =  etJ[i][j]-etM[i]
        if tof_ < 0:
            tof_ = 1*86400
        l2.append(pk.lambert_problem(
            r1 = stateM[i][0:3]*1000,
            r2 = np.asarray(stateJ[i][j][0:3])*1000,
            tof = tof_,
            mu = pk.MU_SUN
        ))

l2 = [l2[i:i+Δ] for i in range(Δ)]
print(l2, type(l2), type(l2[0][0]))

idx = 5
idx2 = 5

l_ = l[idx]
l2_ = l2[idx][idx2]
etM_ = etM[idx]
etJ_ = etJ[idx][idx2]

fig = plt.figure()
ax = fig.gca(projection = '3d')

# PLOTTING BODIES
pl = pk.planet.jpl_lp('earth')
pk.orbit_plots.plot_planet(pl, t0 = pk.epoch(float(spk.et2utc(etE, 'J', 10)[3:]) - 2451544.5), axes = ax, color='b')
pl = pk.planet.jpl_lp('mars')
pk.orbit_plots.plot_planet(pl, t0 = pk.epoch(float(spk.et2utc(etM_, 'J', 10)[3:]) - 2451544.5), axes = ax, color='r')
pl = pk.planet.jpl_lp('jupiter')
pk.orbit_plots.plot_planet(pl, t0 = pk.epoch(float(spk.et2utc(etJ_, 'J', 10)[3:]) - 2451544.5), axes = ax, color='g')

# PLOTTING LAMBERT ARC
pk.orbit_plots.plot_lambert(l_, axes = ax, color='c')
pk.orbit_plots.plot_lambert(l2_, axes = ax, color='c')

# FINDING C3 & Title
v0 = np.array(l_.get_v1()[0])
C3 = np.linalg.norm(stateE[3:6] - v0/1000)**2
tof = l_.get_tof()/84600
plt.title("LAUNCH C3: {:.2f} $km^2/s^2$\nTIME OF FLIGHT: {:.2f} DAYS".format(C3, tof))

ax_a = plt.axes([0.25, 0.08, 0.5, 0.025])
ax_b = plt.axes([0.25, 0.05, 0.5, 0.025])
sl_1 = Slider(ax_b, 'Arc 1 #', 0, len(l)-1, valstep=1, valinit=idx)
sl_2 = Slider(ax_a, 'Arc 2 #', 0, len(l)-1, valstep=1, valinit=0)

def update(val):
    idx = val
    idx2 = val
    ax.cla()
    l_ = l[idx]
    l2_ = l2[idx][idx2]
    etM_ = etM[idx]
    etJ_ = etJ[idx2][idx]

    # PLOTTING BODIES
    pl = pk.planet.jpl_lp('earth')
    pk.orbit_plots.plot_planet(pl, t0 = pk.epoch(float(spk.et2utc(etE, 'J', 10)[3:]) - 2451544.5), axes = ax, color='b')
    pl = pk.planet.jpl_lp('mars')
    pk.orbit_plots.plot_planet(pl, t0 = pk.epoch(float(spk.et2utc(etM_, 'J', 10)[3:]) - 2451544.5), axes = ax, color='r')
    pl = pk.planet.jpl_lp('jupiter')
    pk.orbit_plots.plot_planet(pl, t0 = pk.epoch(float(spk.et2utc(etJ_, 'J', 10)[3:]) - 2451544.5), axes = ax, color='g')

    # PLOTTING LAMBERT ARC
    pk.orbit_plots.plot_lambert(l_, axes = ax, color='c')
    pk.orbit_plots.plot_lambert(l2_, axes = ax, color='c')
    v0 = np.array(l_.get_v1()[0])
    C3 = np.linalg.norm(stateE[3:6] - v0/1000)**2
    tof = l_.get_tof()/84600
    ax.set_title("LAUNCH C3: {:.2f} $km^2/s^2$\nTIME OF FLIGHT: {:.2f} DAYS".format(C3, tof))
    axisEqual3D(ax)

    fig.canvas.draw_idle()

sl_1.on_changed(update)


axisEqual3D(ax)
plt.show()

# # Fixing random state for reproducibility
# np.random.seed(19680801)


# fig, ax = plt.subplots()


# ax.plot(np.random.rand(12), np.random.rand(12), 'go')
# xl = ax.set_xlabel('easy come, easy go')
# ax.set_title('Press a key')
# print(ax.__dict__.keys())
# plt.show()