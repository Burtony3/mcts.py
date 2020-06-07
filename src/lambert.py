import pykep as pk
import numpy as np
import SpicyPy as spk

def lambert(p1, e1, p2, e2):
    # Get Planetary Positions/Velocities
    s1, _ = spk.spkezr(p1, e1, 'J2000', 'NONE', '0')
    s2, _ = spk.spkezr(p2, e2, 'J2000', 'NONE', '0')

    # Time of Flight
    tof = e2 - e1

    # Performing Lambert
    l = pk.lambert_problem(r1  = (s1[0], s1[1], s1[2]),
                           r2  = (s2[0], s2[1], s2[2]),
                           tof = tof,
                           mu  = 1.327e11,
                           max_revs = 1)

    # Calculating V∞'s
    v1 = np.array(l.get_v1()[0])
    v2 = np.array(l.get_v2()[0])
    vInfO = v1 - np.array(s1[3:5])
    vInfI = v2 - np.array(s2[3:5])

    # Returning Variables
    return vInfO , vInfI
