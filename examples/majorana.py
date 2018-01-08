from matplotlib import pyplot
import kwant
import numpy
import tinyarray
import scipy.sparse.linalg

# Python >= 3.3 provides SimpleNamespace in the
# standard library so we can simply import it:
# >>> from types import SimpleNamespace
# (Kwant does not yet support Python 3.)
class SimpleNamespace(object):
    """A simple container for parameters."""
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

s_0 = numpy.identity(2)
s_z = numpy.array([[1, 0], [0, -1]])
s_x = numpy.array([[0, 1], [1, 0]])
s_y = numpy.array([[0, -1j], [1j, 0]])

tau_z = tinyarray.array(numpy.kron(s_z, s_0))
tau_x = tinyarray.array(numpy.kron(s_x, s_0))
sigma_z = tinyarray.array(numpy.kron(s_0, s_z))
tau_zsigma_x = tinyarray.array(numpy.kron(s_z, s_x))

def onsite(site, p):
    return tau_z * (p.mu - 2 * p.t) + \
        sigma_z * p.B + tau_x * p.Delta

def hopping(site0, site1, p):
    return tau_z * p.t + 1j * tau_zsigma_x * p.alpha

def make_system(l=70):
    sys = kwant.Builder()
    lat = kwant.lattice.chain()
    sys[(lat(x) for x in range(l))] = onsite
    sys[lat.neighbors()] = hopping
    return sys.finalized()

sys = make_system()

# Calculate and plot lowest eigenenergies in B-field.
B_values = numpy.linspace(0, 0.6, 80)
energies = []
params = SimpleNamespace(
    t=1, mu=-0.1, alpha=0.05, Delta=0.2)
for params.B in B_values:
    H = sys.hamiltonian_submatrix(
        args=[params], sparse=True)
    H = H.tocsc()
    eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
    energies.append(numpy.sort(eigs[0]))
pyplot.plot(B_values, energies)
pyplot.show()
