r"""
Base class for a Landau-Ginzburg orbifold

AUTHORS:
- Daichi Mukai (2017-01)
"""

from sage.arith.misc import lcm
from sage.matrix.special import diagonal_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.functional import symbolic_sum
from sage.misc.misc_c import prod
from sage.modules.free_module import VectorSpace
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.universal_cyclotomic_field import E

import group_patch

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

    @cached_method
    def central_charge(self):
        r"""
        Return central charge of Landau-Ginzburg orbifold
        """
        return symbolic_sum([ZZ(1) - ZZ(2)*q for q in self.charges()])

    @cached_method
    def __fix__(self, g):
        r"""
        Return fixed subspace w.r.t g
        """
        matrix = g.matrix()
        return (matrix - matrix.parent().one()).kernel()

    @cached_method
    def __block_type__(self):
        return [self.weights().count(w) for w in self.__block_weights__()]

    @cached_method
    def __block_weights__(self):
        weights_no_dup = []
        for w in self.weights():
            if not weights_no_dup.__contains__(w):
                weights_no_dup.append(w)
        return weights_no_dup

    @cached_method
    def __block_subspace__(self):
        V = VectorSpace(CyclotomicField(), len(self.weights()))
        block_type = self.__block_type__()
        S = []
        for i in range(len(block_type)):
            basis = []
            for j in range(block_type[i]):
                offset = symbolic_sum(block_type[0:i])
                base = [0]*len(self.weights())
                base[offset+j] = 1
                basis.append(V(base))
            S.append(V.subspace(basis))
        return S

    @cached_method
    def poincare_series(self):
        conj = self.group().conjugacy_class_representatives()
        return symbolic_sum([self.poincare_series_sector(g) for g in conj])

    def poincare_series_sector(self, g):
        fix = self.__fix__(g)
        centralizer = None if fix.rank() == 0 else self.group().centralizer(g)

        left, right = self.__vacuum_shift__(g)
        l = PolynomialRing(ZZ, 'u, v')
        u, v = l.gens()
        return self.__poincare_series_sector_sum__(fix, centralizer) * u**left * v**right

    @cached_method
    def __vacuum_shift__(self, g):
        def _log_2pii_(x):
            i = 0
            order = x.multiplicative_order()
            z = E(order)
            for i in range(order):
                if z**i == x:
                    return ZZ(i)/order
            raise ArithmeticError
        S = self.__block_subspace__()
        ev_g = group_patch.eigenvalues(g.matrix())
        plus = ZZ(0)
        minus = ZZ(0)
        for i in range(len(S)):
            evs = group_patch.eigenvalues(g.matrix().restrict(S[i]))
            plus += sum([self.__block_weights__()[i]-self.degree()/2 for j in range(len(evs)) if evs[j] == 1])
            minus += sum([_log_2pii_(evs[j])*self.degree()-self.degree()/2 for j in range(len(evs)) if evs[j] != 1])
        plus += self.central_charge()*self.degree()/2
        return (plus+minus, plus-minus)

    @cached_method
    def __poincare_series_sector_sum__(self, fix, centralizer):
        l = PolynomialRing(ZZ, 'u, v')
        u, v = l.gens()
        ret = l.zero()
        if fix.rank() == 0:
            assert centralizer == None
            ret = l.one()
        else:
            k_poly = PolynomialRing(CyclotomicField(), 't')
            t = k_poly.gens()[0]
            tmp_ps = k_poly.zero()
            S = self.__block_subspace__()
            _coerce_ = CyclotomicField()._coerce_
            for h in centralizer:
                h_decomposition = [h.matrix().restrict(fix.intersection(s)) for s in S]
                tmp = k_poly.one()
                for i in range(len(h_decomposition)):
                    h_sub = h_decomposition[i]
                    if h_sub.rank() == 0:
                        next
                    evs = [_coerce_(e) for e in sorted(group_patch.eigenvalues(h_sub))]
                    tmp *= prod([(ev-(t)**(self.degree()-self.__block_weights__()[i]))/(1-ev*(t)**(self.__block_weights__()[i])) for ev in evs])
                tmp_ps += tmp
            ret = tmp_ps(t=u*v)/len(centralizer)
        return ret
