r"""
Patches for Group operations

AUTHORS:
- Daichi Mukai (2017-01)
"""

from sage.groups.matrix_gps.finitely_generated  import FinitelyGeneratedMatrixGroup_gap
from sage.misc.cachefunc import cached_method

@cached_method
def __centralizer__(self, element):
    g = element.gap()
    C = self.gap().Centralizer(g)
    return self.subgroup(C.GeneratorsOfGroup())
FinitelyGeneratedMatrixGroup_gap.centralizer = __centralizer__
