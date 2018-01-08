import matplotlib.pyplot
import kwant

lat = kwant.lattice.square()
sys = kwant.Builder()

# Add sites.
sys[lat(0, 0)] = 1.5
sys[lat(1, 0)] = 1.5
sys[lat(0, 1)] = 1.5

# Add hoppings.
sys[lat(0, 0), lat(1, 0)] = 2j
sys[lat(0, 1), lat(1, 0)] = 2j

kwant.plot(sys)
