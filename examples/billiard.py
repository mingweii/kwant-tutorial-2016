import kwant

# Build the scattering region.
sys = kwant.Builder()
sqlat = kwant.lattice.square()

def stadium(position):
    x, y = position
    x = max(abs(x) - 70, 0)
    return x**2 + y**2 < 100**2

sys[sqlat.shape(stadium, (0, 0))] = 4
sys[sqlat.neighbors()] = -1

# Build and attach the leads, finalize system.
lead_symmetry = kwant.TranslationalSymmetry([0, -1])
for start, end in [(-90, -60), (0, 30)]:
    lead = kwant.Builder(lead_symmetry)
    lead[(sqlat(x, 0) for x in range(start, end))] = 4
    lead[sqlat.neighbors()] = -1
    sys.attach_lead(lead)

sys = sys.finalized()

# Compute and plot observables.
energies = [0.5 + 1e-4 * i for i in range(300)]
conductances = [kwant.smatrix(sys, en).transmission(1, 0)
                for en in energies]

local_dos = kwant.ldos(sys, energy=.2)

from matplotlib import pyplot
pyplot.plot(energies, conductances)
pyplot.show()
kwant.plotter.map(sys, local_dos, num_lead_cells=10)
