import tinyarray as ta
import math
from matplotlib import pyplot
from matplotlib.patches import Circle
import kwant

def circle(pos):
    x, y = pos
    return x**2 + y**2 < 20

lat = kwant.lattice.triangular()
sys = kwant.Builder()
sys[lat.shape(circle, (0, 0))] = 0
sys[lat.neighbors()] = 1

lead_dirs = [lat.vec((-3, 1)), lat.vec((0, -1)), lat.vec((2, -1))]
for d in lead_dirs:
    lead = kwant.Builder(kwant.TranslationalSymmetry(d))
    lead[lat.wire((0, 0), 2.1)] = 0
    lead[lat.neighbors()] = 1
    sys.attach_lead(lead)

fig = kwant.plot(sys, show=False)
ax = pyplot.gca()
pyplot.text(-2, 5.5, 'scattering region', size=15)
pyplot.text(-10, -1, 'lead 0', color='red', size=15)
pyplot.text(-3, -7.7, 'lead 1', color='red', size=15)
pyplot.text(5.5, 0, 'lead 2', color='red', size=15)

for dir, offset in zip(lead_dirs, [10.5, 8, 8.6]):
    dir = dir / math.sqrt(ta.dot(dir, dir))
    for i in [0, 0.4, 0.8]:
        ax.add_artist(Circle(dir * (offset + i), 0.06, fc='k'))

pyplot.axis('off')
pyplot.xlim((-11, 9))
pyplot.ylim((-8, 6))
fig.tight_layout()
fig.set_size_inches(*(5, 3.75))
fig.savefig("tbsys.pdf")
