import math
from cmath import exp
import numpy
from matplotlib import pyplot

import kwant
from kwant.digest import gauss

def hopping(sitei, sitej, phi, salt):
    xi, yi = sitei.pos
    xj, yj = sitej.pos
    return -exp(-0.5j * phi * (xi - xj) * (yi + yj))

def onsite(site, phi, salt):
    return 0.05 * gauss(repr(site), salt) + 4

def make_system(L=50):
    def central_region(pos):
        x, y = pos
        return -L < x < L and \
            abs(y) < L - 37.5 * math.exp(-x**2 / 12**2)

    lat = kwant.lattice.square()
    sys = kwant.Builder()

    sys[lat.shape(central_region, (0, 0))] = onsite
    sys[lat.neighbors()] = hopping

    sym = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym)
    lead[(lat(0, y) for y in range(-L + 1, L))] = 4
    lead[lat.neighbors()] = hopping

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

sys = make_system()
energy = 0.15

# Calculate and plot QHE conductance plateaus.
reciprocal_phis = numpy.linspace(4, 50, 200)
conductances = []
for phi in 1 / reciprocal_phis:
    smatrix = kwant.smatrix(sys, energy, args=[phi, ""])
    conductances.append(smatrix.transmission(1, 0))
pyplot.plot(reciprocal_phis, conductances)
pyplot.show()

# Calculate and plot a QHE edge state.
def density(sys, energy, args, lead_nr):
    wf = kwant.wave_function(sys, energy, args)
    return (abs(wf(lead_nr))**2).sum(axis=0)

d = density(sys, energy, [1/40.0, ""], 0)
kwant.plotter.map(sys, d)
