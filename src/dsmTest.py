import spiceypy as spk
import numpy as np
from matplotlib import pyplot as plt
import math
from numpy.linalg import norm

spk.furnsh("../data/spk/naif0009.tls")
spk.furnsh("../data/spk/de438.bsp")
spk.furnsh("../data/spk/gm_de431.tpc")
spk.furnsh("../data/spk/pck00010.tpc")

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

AU = 149597870700/1000

utc = "Jan 01, 2020"

# tof = 365.2*2 + (49/360)*365.2
tof = 715.537

et0 = spk.utc2et(utc)

et1 = et0 + tof*86400

pos0 = spk.spkezr('3', et0, "ECLIPJ2000", "NONE", '0')[0][:3]/AU

pos1 = spk.spkezr('3', et1, "ECLIPJ2000", "NONE", '0')[0][:3]/AU

theta = math.acos(np.dot(pos0, pos1)/(norm(pos0)*norm(pos1))) * 180/math.pi

print("Angle between vectors: {} Degrees\n   Should be: 49 Degrees".format(theta))

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.view_init(azim=0, elev=90)
tau = 365.2*86400
et1 = et0 + tau
et = np.linspace(et0, et1, 250)
x = []
y = []
z = []
for et_ in et:
    r = spk.spkezr('3', et_, "ECLIPJ2000", "NONE", '0')[0][:3]/AU
    x.append(r[0])
    y.append(r[1])
    z.append(r[2])
ax.plot3D(x, y, z, 'gray')
ax.scatter(pos0[0], pos0[1], pos0[2], label="Launch")
ax.scatter(pos1[0], pos1[1], pos1[2], label="Flyby")
ax.legend()
# axisEqual3D(ax)
plt.show()