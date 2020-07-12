import pykep as pk
import numpy as np
import SpicyPy as spk

def lambert(p0, e0, p1, e1):
    # Get Planetary Positions/Velocities
    s0, _ = spk.spkezr(p0, e0, 'J2000', 'NONE', '0')
    s1, _ = spk.spkezr(p1, e1, 'J2000', 'NONE', '0')

    # Time of Flight
    tof = e1 - e0

    # Performing Lambert
    l = pk.lambert_problem(r0  = (s0[0], s0[1], s0[2]),
                           r1  = (s1[0], s1[1], s1[2]),
                           tof = tof,
                           mu  = 1.327e01,
                           max_revs = 1)

    # Calculating V∞'s
    v0 = np.array(l.get_v0()[0])
    v1 = np.array(l.get_v1()[0])
    vInfO = v0 - np.array(s0[3:5])
    vInfI = v1 - np.array(s1[3:5])

    # Returning Variables
    return vInfO , vInfI
