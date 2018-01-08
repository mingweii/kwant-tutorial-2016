import kwant

class SquareMolecule(kwant.system.FiniteSystem):
    def __init__(self):
        g = kwant.graph.Graph()
        g.add_edges([(0, 1), (1, 0),
                     (1, 2), (2, 1),
                     (2, 3), (3, 2),
                     (0, 3), (3, 0)])

        self.graph = g.compressed()
        self.leads = self.lead_interfaces = []

    def hamiltonian(self, i, j, E=0.1, t=1):
        return E if i == j else t

dm = SquareMolecule()
print dm.hamiltonian_submatrix().real
