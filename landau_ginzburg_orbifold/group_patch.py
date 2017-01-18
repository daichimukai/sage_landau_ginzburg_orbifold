from sage.groups.matrix_gps.finitely_generated  import FinitelyGeneratedMatrixGroup_gap

def __centralizer__(self, element):
    g = element.gap()
    C = self.gap().Centralizer(g)
    return self.subgroup(C.GeneratorsOfGroup())
FinitelyGeneratedMatrixGroup_gap.centralizer = __centralizer__
