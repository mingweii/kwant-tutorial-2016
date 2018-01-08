from __future__ import division

from math import sqrt
import tinyarray as ta
import numpy as np
from scipy.spatial import cKDTree as KDTree
# from itertools import izip <<--- python 2 only
from concurrent import futures  # Available from futures package.
import kwant


# Pauli matrices
zero = ta.array(((0, 0), (0, 0)), complex)
s0 = ta.array(((1, 0), (0, 1)), complex)
sx = ta.array(((0, 1), (1, 0)), complex)
sy = ta.array(((0, -1j), (1j, 0)), complex)
sz = ta.array(((1, 0), (0, -1)), complex)


def make_system(onsite_sys, hop_sys, L, W, leads, onsite_g, hop_g):
    """Construct the a lateral spin valve.

    Parameters:
    -----------
    onsite_sys : Onsite Hamiltonian in the system.
    hop_sys : Hopping in the system.
    L : Half the system length.
    W : Half the system width.
    leads : Sequence of 5-tuples specifying the top leads.
            Each tuple has the form (x_start, x_end, onsite, hop, vhop).
    onsite_g : Onsite Hamiltonian in the in-plane graphene leads
    hop_g : Hopping the the in-plane graphene leads
    """
    # Define lattices
    # Honeycomb lattice (only z=0 sites are used.)
    honeycomb = kwant.lattice.general([(1, 0, 0), (0.5, 0.5 * sqrt(3), 0), (0, 0, 1)],
                                [(0, 0, 0), (0, 1 / sqrt(3), 0)])
    # Cubic lattice
    cubic = kwant.lattice.general(np.identity(3))

    def rect(x0, y0, x1, y1, z=None):
        """Return a tuple (f, origin) of arguments for Kwant's shape function.

        The shape is a rectangle in the `xy`-plane.  If `z` differs from
        `None`, that coordinate will be restricted to the given value.
        """
        def inside(pos):
            px, py, pz = pos
            return (z is None or pz == z) and x0 <= px <= x1 and y0 <= py <= y1

        return inside, (0.5 * (x0 + x1), 0.5 * (y0 + y1), 0 if z is None else z)

    # Build system.
    sys = kwant.Builder()
    sys[honeycomb.shape(*rect(-L, -W, L, W, 0))] = onsite_sys
    graphene_sites = list(sys.sites())
    kdt = KDTree(np.array([site.pos for site in graphene_sites])[:, :2])

    # Add leads.
    def add_lead(x0, y0, x1, y1, onsite_lead, hop_lead, hop_inter):
        """Create and attach a vertical lead to the system.

        One unit cell of the lead is added to the system so that the hopping
        between the system and the lead can be controlled.
        """
        # Add one unit cell of the (yet to be constructed) lead and connect it
        # to the honeycomb lattice of the system.
        for a in sys.expand(cubic.shape(*rect(x0, y0, x1, y1, 1))):
            sys[a] = onsite_lead
            x, y, z = a.pos
            dist, idx = kdt.query((x, y))
            sys[a, graphene_sites[idx]] = hop_inter
            sys[cubic.neighbors()] = hop_lead

        # Construct the lead and attach it to the system.
        sym = kwant.TranslationalSymmetry([0, 0, 1])
        lead = kwant.Builder(sym)
        lead[cubic.shape(*rect(x0, y0, x1, y1))] = onsite_lead
        lead[cubic.neighbors()] = hop_lead
        sys.attach_lead(lead)

    def add_graphene_lead(x0, y0, x1, y1, onsite, hop, direction):
        """Add a graphene lead."""
        sym = kwant.TranslationalSymmetry([direction, 0, 0])
        lead = kwant.Builder(sym)
        lead[honeycomb.shape(*rect(x0, y0, x1, y1,0))] = onsite
        lead[honeycomb.neighbors()] = hop
        sys.attach_lead(lead)

    for x0, x1, onsite_lead, hop_lead, hop_inter in leads:
        add_lead(x0, -W, x1, W, onsite_lead, hop_lead, hop_inter)
    add_graphene_lead(-W, -W, W, W, onsite_g, hop_g, +1)
    add_graphene_lead(-W, -W, W, W, onsite_g, hop_g, -1)

    # Set the in-plane hoppings of the graphene sheet.
    sys[honeycomb.neighbors()] = hop_sys

    return sys


def conductance_matrix(sys, energy, args=()):
    n = len(sys.leads)
    s = kwant.smatrix(sys, energy, args)
    cond = np.array([[s.transmission(i, j) for j in range(n)] for i in
                     range(n)])
    # Correct the reflection blocks, so that rows and columns sum to zero.
    cond -= np.diag(cond.sum(axis=0))
    return cond


class SimpleNamespace(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def onsite(s, p):
    return (kwant.digest.uniform(s.pos, 'E') - 0.5) * p.dis * s0 + p.H * sx


def hop_lead_left_interface(a,b, p):
    return p.tsL / 2 * s0 + p.betaL * p.tsL / 2 * sz


def hop_lead_right_interface(a,b, p):
    return p.tsR / 2 * s0 + p.betaR * p.tsR / 2 * sz


P = 4
sys = make_system(onsite, s0, 6*P, 8*P,
                  [(-2*P, -P, zero, s0, hop_lead_left_interface),
                   (P, 2*P, zero, s0, hop_lead_right_interface)],
                  zero, s0).finalized()


def calc(H):
    #### Non-local conductance ####
    # We output the voltage difference in the right electrodes as a function of
    # current between the left electrodes.
    args = SimpleNamespace(dis=0.4, tsL=0.2, betaL=0.5, tsR=0.2, betaR=0.5)
    energy = 0.75

    args.H = H
    # In order to calculate the nonlocal conductance we eliminate one row and
    # one column from the condutcance matrix. This amounts to setting the
    # corresponding voltage to zero and using current conservation to calculate
    # the current through the last terminal.
    cm_P = conductance_matrix(sys, energy, [args])[:-1, :-1]
    # We then set the current to be 1 in the lead 0, -1 in lead 2, and
    # calculate the voltage in lead 1 (so V_1 - V_3 since V_3 = 0).
    nlsP = np.linalg.solve(cm_P, [1, 0, -1])[1]

    args.betaR = -args.betaR
    cm_AP = conductance_matrix(sys, energy, [args])[:-1, :-1]
    nlsAP = np.linalg.solve(cm_AP, [1, 0, -1])[1]

    return nlsP, nlsAP


def main():
    H_values = np.linspace(-0.2, 0.2, 192)
    with futures.ProcessPoolExecutor() as executor:
        for H, result in zip(H_values, executor.map(calc, H_values)):
            print ( H, result[0], result[1])


if __name__ == '__main__':
    main()
