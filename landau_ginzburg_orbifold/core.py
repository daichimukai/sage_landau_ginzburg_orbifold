r"""
Base class for a Landau-Ginzburg orbifold

AUTHORS:
- Daichi Mukai (2017-01)
"""

from sage.arith.misc import lcm
from sage.matrix.special import diagonal_matrix
from sage.misc.cachefunc import cached_method
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.universal_cyclotomic_field import E

class LandauGinzburgOrbifold:
    def __init__(self, group, charges):
        if group == None or charges == None:
            raise ValueError

        self._group_ = group
        self._charges_ = charges

    def group(self):
        return self._group_

    def charges(self):
        r"""
        Return charges of polynomial
        """
        return self._charges_

    @cached_method
    def degree(self):
        r"""
        Return degree of polynomial
        """
        return lcm([q.denominator() for q in self.charges()])

    @cached_method
    def weights(self):
        r"""
        Return list of charges of variables
        """
        d = self.degree()
        return [q*d for q in self.charges()]

    @cached_method
    def exponential_grading_operator(self):
        r"""
        Return the exponential grading operator
        """
        d = self.degree()
        return diagonal_matrix(CyclotomicField(), [E(d)**w for w in self.weights()])
