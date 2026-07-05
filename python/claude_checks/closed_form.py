"""Closed-form candidate for the single-column (i=0) non-free deviation:

constants: c_m(w) = binom(m - w, 2) mod 2 for m >= 3, else 0
rule:      F_n(u)(V) = sum_{A subset V} c_{n-|A|}(u - |A cap [0,u)|)
                     = sum_{j<=p, l<=q} C(p,j) C(q,l) c_{n-j-l}(u-j),
           where p = |V cap [0,u)|, q = |V cap (u,n]|  (type-only!)
family:    X[n] = sum_u sum_{partitions (V,W) of c(u)} F_n(u)(V) . V ot W

Check: symmetry, single-column chain identity bd X[n] = sum_j delta_j X[n-1],
irreducibility (by construction), for n <= 8. Then full construction checks.
"""
from math import comb
from itertools import combinations
from cupi import bd, T, one_plus_T, cofaces_sum, check_identity, canonical, is_free, is_irreducible

def c_const(m, w):
    if m < 3 or w < 0 or w > m:
        return 0
    return comb(m - w, 2) % 2

def f_type(n, u, p, q):
    s = 0
    for j in range(p + 1):
        for l in range(q + 1):
            if comb(p, j) % 2 and comb(q, l) % 2:
                s ^= c_const(n - j - l, u - j)
    return s

def X_level(n):
    out = set()
    for u in range(n + 1):
        cu = [x for x in range(n + 1) if x != u]
        below = [x for x in cu if x < u]
        for r in range(len(cu) + 1):
            for V in combinations(cu, r):
                Vs = frozenset(V)
                p = len([x for x in V if x < u])
                q = r - p
                if f_type(n, u, p, q):
                    out ^= {(Vs, frozenset(cu) - Vs)}
    return out

NMAX = 8
X = {n: X_level(n) for n in range(NMAX + 1)}

print("level sizes:", {n: len(X[n]) for n in range(NMAX + 1)})
ok_sym = all(one_plus_T(X[n]) == set() for n in range(NMAX + 1))
print("symmetric at all levels:", ok_sym)
ok_irr = all(is_irreducible(X[n]) for n in range(NMAX + 1))
print("irreducible at all levels:", ok_irr)

ok_chain = True
for n in range(1, NMAX + 1):
    defect = bd(X[n], n) ^ cofaces_sum(X[n - 1], n)
    if defect:
        ok_chain = False
        print(f"chain identity FAILS at n={n}: defect size {len(defect)}")
print(f"chain identity bd X[n] = sum delta X[n-1] for n <= {NMAX}:", ok_chain)

# full construction: canonical + X at index 0
def D(i, n):
    base = canonical(i, n)
    if i == 0 and n <= NMAX:
        return base ^ X[n]
    return base

bad = []
for n in range(7):
    for i in range(8):
        d = check_identity(D, i, n)
        if d:
            bad.append((i, n))
print("full construction identities n <= 6:", "all pass" if not bad else f"FAIL {bad}")
print("non-degenerate (Delta_0[0] = empty ot empty):", D(0, 0) == canonical(0, 0))
nf = [(0, n) for n in range(NMAX + 1) if not is_free(D(0, n), 0, n)]
print("freeness fails at (i=0, n):", [n for _, n in nf])
print("first deviation level:", min(n for n in range(NMAX + 1) if X[n]))
