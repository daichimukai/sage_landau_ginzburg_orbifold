r"""
Patches for Group operations

AUTHORS:
- Daichi Mukai (2017-01)
"""

from sage.arith.misc import lcm
from sage.groups.matrix_gps.finitely_generated  import FinitelyGeneratedMatrixGroup_gap
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import CyclotomicField

@cached_method
def __centralizer__(self, element):
    g = element.gap()
    C = self.gap().Centralizer(g)
    return self.subgroup(C.GeneratorsOfGroup())
FinitelyGeneratedMatrixGroup_gap.centralizer = __centralizer__

def eigenvalues(element):
    def __multiplicative_order__(element):
        g = element
        one = element.parent().one()
        i = 1
        while True:
            if g == one:
                return ZZ(i)
            g = g * element
            i += 1

    order = lcm([value.conductor() for value in element.dict().values()] + [__multiplicative_order__(element)])
    field = CyclotomicField(order)
    return element.change_ring(field).eigenvalues()
