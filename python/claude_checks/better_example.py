"""Verify the 'better' branch non-free example:
Delta_0[2] += zeta, zeta = (1+T)(empty ot {02} + empty ot {01} + {2} ot {0} + {1} ot {0})"""
from cupi import (canonical, T, one_plus_T, bd, check_identity,
                  is_irreducible, is_free, cofaces_sum)

f = frozenset
half = {(f(), f({0, 2})), (f(), f({0, 1})), (f({2}), f({0})), (f({1}), f({0}))}
zeta = half ^ T(half)

print(f"zeta is a cycle: {bd(zeta, 2) == set()}")
print(f"zeta symmetric: {one_plus_T(zeta) == set()}")

def better(i, n):
    if (i, n) == (0, 2):
        return canonical(0, 2) ^ zeta
    return canonical(i, n)

bad = []
for n in range(8):
    for i in range(9):
        d = check_identity(better, i, n)
        if d:
            bad.append((i, n, len(d)))
print(f"failing identities up to n=7: {bad if bad else 'NONE - valid construction'}")
t = better(0, 2)
print(f"Delta_0[2]: irreducible={is_irreducible(t)}, free={is_free(t, 0, 2)}")
print(f"Delta_0[0] = canonical (non-degenerate): {better(0,0) == canonical(0,0)}")

# also check simpliciality of the perturbed component
from simplicial_check import codegeneracy_subset
simp = True
for j in range(2):
    out = set()
    for V, W in t:
        V2, W2 = codegeneracy_subset(V, j), codegeneracy_subset(W, j)
        if V2 is not None and W2 is not None:
            out ^= {(V2, W2)}
    if out:
        simp = False
print(f"perturbed Delta_0[2] simplicial: {simp}")

# key structural fact: why does this work while bd zeta_{1} fails?
push = cofaces_sum(zeta, 3)
print(f"sum_j P(delta_j)^ot2 (zeta) = 0: {push == set()}  <-- the missing condition")
