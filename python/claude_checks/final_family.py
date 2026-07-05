"""THE closed-form counterexample family (single column, i=0):

constants: c_m(w) = binom(m,2) + binom(w,2) mod 2 for m >= 3, else 0
           (equivalently c = a + Ta for the seed a = xy + y^2)
rules:     F_n(u)(V) = sum_{j<=p, l<=q} C(p,j) C(q,l) c_{n-j-l}(u-j),
           p = |V cap [0,u)|, q = |V cap (u,n]|
family:    X[n] = sum_u sum_{(V,W) partition of c(u)} F_n(u)(V) V ot W

Full verification through n = 10."""
from math import comb
from itertools import combinations
from cupi import bd, one_plus_T, cofaces_sum, check_identity, canonical, is_free, is_irreducible

def c_const(m, w):
    if m < 3 or w < 0 or w > m:
        return 0
    return (comb(m, 2) + comb(w, 2)) % 2

def f_type(n, u, p, q):
    s = 0
    for j in range(p + 1):
        if comb(p, j) % 2 == 0: continue
        for l in range(q + 1):
            if comb(q, l) % 2 == 0: continue
            s ^= c_const(n - j - l, u - j)
    return s

def X_level(n):
    out = set()
    for u in range(n + 1):
        cu = [x for x in range(n + 1) if x != u]
        for r in range(len(cu) + 1):
            for V in combinations(cu, r):
                p = len([x for x in V if x < u])
                if f_type(n, u, p, r - p):
                    Vs = frozenset(V)
                    out ^= {(Vs, frozenset(cu) - Vs)}
    return out

NMAX = 10
X = {n: X_level(n) for n in range(NMAX + 1)}
print("level sizes:", {n: len(X[n]) for n in range(NMAX + 1)})
print("first deviation level:", min(n for n in range(NMAX + 1) if X[n]))
print("X[3] check (= zeta_c(0) + zeta_c(1)): ", end="")
from cupi import partition_chain
target = partition_chain(frozenset({1,2,3})) ^ partition_chain(frozenset({0,2,3}))
print(X[3] == target)

ok_sym = all(one_plus_T(X[n]) == set() for n in range(NMAX + 1))
ok_irr = all(is_irreducible(X[n]) for n in range(NMAX + 1))
ok_chain = all(bd(X[n], n) ^ cofaces_sum(X[n - 1], n) == set() for n in range(1, NMAX + 1))
print(f"symmetric: {ok_sym}, irreducible: {ok_irr}, chain identity n<=({NMAX}): {ok_chain}")

def D(i, n):
    base = canonical(i, n)
    if i == 0 and n <= NMAX:
        return base ^ X[n]
    return base

bad = [(i, n) for n in range(7) for i in range(8) if check_identity(D, i, n)]
print("full construction (canonical + X) identities n <= 6:", "all pass" if not bad else f"FAIL {bad}")
print("non-degenerate:", D(0, 0) == canonical(0, 0))
print("NOT free at (0,n) for n in:", [n for n in range(NMAX + 1) if not is_free(D(0, n), 0, n)])
