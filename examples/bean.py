import matplotlib.pyplot
import kwant

lat = kwant.lattice.honeycomb()
sys = kwant.Builder()

def bean(pos):
    x, y = pos
    rr = x**2 + y**2
    return rr**2 < 15 * y * rr + x**2 * y**2

sys[lat.shape(bean, (0, 1))] = 0.5
sys[lat.neighbors(1)] = 1

kwant.plot(sys)
