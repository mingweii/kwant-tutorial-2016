import matplotlib.pyplot
import kwant

lat = kwant.lattice.square()
sys = kwant.Builder()

def ring(pos):
    x, y = pos
    return 7**2 <= x**2 + y**2 < 13**2

sys[lat.shape(ring, (10, 0))] = 0
sys[lat.neighbors()] = 1

sym = kwant.TranslationalSymmetry((-2, 1))
lead = kwant.Builder(sym)
lead[(lat(x, 0) for x in range(-6, 6))] = 0
lead[lat.neighbors()] = 1
sys.attach_lead(lead)
sys.attach_lead(lead, lat(0, 0))

kwant.plot(sys)
