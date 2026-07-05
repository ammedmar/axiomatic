"""
Verify cup-i construction identities in the P-model of sec/reformulation.tex.

A chain in P(simplex^n)^{ot 2} is a set of pairs (V, W) with V, W subsets of
{0,...,n} (subsets of removed face indices), coefficients in F_2.

Identity to check at (i, n):
    bd Delta_i[n] + sum_j P(delta_j)^{ot2} Delta_i[n-1] = (1+T) Delta_{i-1}[n]
"""

from itertools import combinations


def add(a, b):
    return a ^ b  # symmetric difference = mod 2 sum of sets of basis elements


def bd_subset(U, n):
    """Boundary of a single subset U in P(simplex^n): sum of U + {ubar}.

    U represents the face d_U[n] of degree n - |U|; if |U| = n (a vertex),
    the boundary is zero.
    """
    if len(U) >= n:
        return set()
    return {frozenset(U | {u}) for u in range(n + 1) if u not in U}


def bd(chain, n):
    """Boundary on tensor product chain."""
    out = set()
    for V, W in chain:
        for V2 in bd_subset(V, n):
            out ^= {(V2, W)}
        for W2 in bd_subset(W, n):
            out ^= {(V, W2)}
    return out


def T(chain):
    return {(W, V) for V, W in chain}


def one_plus_T(chain):
    return chain ^ T(chain)


def coface_subset(U, j):
    """P(delta_j): insert j, shift elements >= j up by one."""
    return frozenset({u if u < j else u + 1 for u in U} | {j})


def cofaces_sum(chain, n_target):
    """sum_{j=0}^{n_target} P(delta_j)^{ot2} applied to chain over simplex^{n_target-1}."""
    out = set()
    for j in range(n_target + 1):
        for V, W in chain:
            out ^= {(coface_subset(V, j), coface_subset(W, j))}
    return out


def canonical(i, n):
    """Canonical cup-i construction in P-model: sum over U in P^n_{n-i} of U^0 ot U^1."""
    if i > n or i < 0:
        return set()
    out = set()
    for U in combinations(range(n + 1), n - i):
        U0 = frozenset(u for k, u in enumerate(U, start=1) if (u + k) % 2 == 0)
        U1 = frozenset(u for k, u in enumerate(U, start=1) if (u + k) % 2 == 1)
        out ^= {(U0, U1)}
    return out


def check_identity(D, i, n):
    """D: function (i, n) -> chain. Returns defect (empty set means identity holds)."""
    lhs = bd(D(i, n), n) ^ cofaces_sum(D(i, n - 1), n) if n > 0 else bd(D(i, n), n)
    rhs = one_plus_T(D(i - 1, n)) if i > 0 else set()
    return lhs ^ rhs


def is_irreducible(chain):
    return all(not (V & W) for V, W in chain)


def is_free(chain, i, n):
    if i == n:
        return True
    return all((W, V) not in chain or (V, W) == (W, V) for V, W in chain)


def partition_chain(U):
    """zeta_U: all (V, W) with V | W = U, V & W = empty."""
    U = list(U)
    out = set()
    for r in range(len(U) + 1):
        for V in combinations(U, r):
            V = frozenset(V)
            out.add((V, frozenset(U) - V))
    return out


NMAX = 5

print("=== 1. Sanity: canonical construction satisfies all identities ===")
ok = True
for n in range(NMAX + 1):
    for i in range(NMAX + 1):
        d = check_identity(canonical, i, n)
        if d:
            ok = False
            print(f"  canonical FAILS at (i={i}, n={n}): defect {d}")
print("  all pass" if ok else "  FAILURE")

print()
print("=== 2. Non-free example: tilde_Delta_0[2] = canonical_0[2] + bd zeta_{1} ===")
zeta1 = partition_chain({1})  # in P(simplex^2)^{ot2}
bd_zeta1 = bd(zeta1, 2)


def tilde(i, n):
    if (i, n) == (0, 2):
        return canonical(0, 2) ^ bd_zeta1
    return canonical(i, n)


t02 = tilde(0, 2)
print(f"  tilde_Delta_0[2] irreducible: {is_irreducible(t02)}")
print(f"  tilde_Delta_0[2] free: {is_free(t02, 0, 2)}")
for n in range(NMAX + 1):
    for i in range(NMAX + 1):
        d = check_identity(tilde, i, n)
        if d:
            print(f"  identity FAILS at (i={i}, n={n}), defect has {len(d)} terms:")
            for V, W in sorted(d, key=lambda p: (sorted(p[0]), sorted(p[1]))):
                print(f"    {tuple(sorted(V))} ot {tuple(sorted(W))}")

print()
print("=== 3. Naive repair: also perturb Delta_0[3] by sum_j P(delta_j)^{ot2} zeta_{1} ===")
push_zeta = cofaces_sum(zeta1, 3)


def tilde2(i, n):
    if (i, n) == (0, 2):
        return canonical(0, 2) ^ bd_zeta1
    if (i, n) == (0, 3):
        return canonical(0, 3) ^ push_zeta
    return canonical(i, n)


bad = []
for n in range(NMAX + 1):
    for i in range(NMAX + 1):
        d = check_identity(tilde2, i, n)
        if d:
            bad.append((i, n, len(d)))
print(f"  failing identities: {bad if bad else 'none up to n=%d' % NMAX}")
t03 = tilde2(0, 3)
print(f"  repaired Delta_0[3] irreducible: {is_irreducible(t03)}")
print(f"  repaired Delta_0[3] free: {is_free(t03, 0, 3)}")

print()
print("=== 4. Transposing-one-element example (i0=1, n0=3): free, not irreducible ===")
i0, n0 = 1, 3


def transposed(i, n):
    if (i, n) == (i0, n0):
        return T(canonical(i0, n0))
    if (i, n) == (i0 - 1, n0 + 1):
        return canonical(i0 - 1, n0 + 1) ^ cofaces_sum(canonical(i0, n0), n0 + 1)
    if (i, n) == (i0 - 1, n0):
        return canonical(i0 - 1, n0) ^ bd(canonical(i0, n0), n0)
    return canonical(i, n)


bad = []
for n in range(NMAX + 2):
    for i in range(NMAX + 2):
        d = check_identity(transposed, i, n)
        if d:
            bad.append((i, n, len(d)))
print(f"  failing identities: {bad if bad else 'none up to n=%d' % (NMAX + 1)}")
for (ii, nn) in [(i0, n0), (i0 - 1, n0 + 1), (i0 - 1, n0)]:
    ch = transposed(ii, nn)
    print(f"  Delta_{ii}[{nn}]: irreducible={is_irreducible(ch)}, free={is_free(ch, ii, nn)}")

print()
print("=== 5. Does ANY irreducible correction c at (0,3) repair the example? ===")
# Need c in P(simplex^3)^{ot2}, degree 3, with:
#   (a) support on irreducible basis elements (V & W = empty),
#   (b) bd c = D  (the 32-term defect),
#   (c) (1+T) c = 0  (so the (1,3) identity is untouched).
# Solve the F_2 linear system by Gaussian elimination over orbit variables.

D = check_identity(tilde, 0, 3)  # the defect
assert one_plus_T(D) == set() or T(D) == D  # D is symmetric

# irreducible degree-3 basis elements of P(simplex^3)^{ot2}: |V|+|W|=3, disjoint
irr_basis = []
for r in range(4):
    for V in combinations(range(4), r):
        rest = [u for u in range(4) if u not in V]
        for W in combinations(rest, 3 - r):
            irr_basis.append((frozenset(V), frozenset(W)))
print(f"  irreducible degree-3 basis size: {len(irr_basis)}")

# symmetric orbit generators: (V,W)+(W,V)
orbits = []
seen = set()
for b in irr_basis:
    if b in seen:
        continue
    tb = (b[1], b[0])
    seen.add(b)
    seen.add(tb)
    orbits.append({b, tb} if tb != b else {b})

# each orbit variable maps to bd(orbit) in degree-2 target space
cols = [bd(set().union(o), 3) for o in orbits]

# Gaussian elimination over F_2: represent chains as frozensets of basis pairs
target = set(D)
pivots = []  # list of (chain, varset)
varsets = [ {k} for k in range(len(orbits)) ]
rows = [ (set(c), {k}) for k, c in enumerate(cols) ]

# reduce columns to echelon form
basis_echelon = []
for chain, vs in rows:
    chain = set(chain)
    vs = set(vs)
    for e_chain, e_vs in basis_echelon:
        pivot = next(iter(sorted(e_chain, key=repr)))
        if pivot in chain:
            chain ^= e_chain
            vs ^= e_vs
    if chain:
        basis_echelon.append((chain, vs))

# now reduce target
t = set(target)
tvs = set()
for e_chain, e_vs in basis_echelon:
    pivot = next(iter(sorted(e_chain, key=repr)))
    if pivot in t:
        t ^= e_chain
        tvs ^= e_vs
if t:
    print("  NO irreducible symmetric c with bd c = D exists (system inconsistent).")
else:
    sol = set()
    for k in tvs:
        sol ^= set().union(orbits[k])
    print(f"  solution found with {len(sol)} terms; verifying...")
    assert bd(sol, 3) == D
    assert one_plus_T(sol) == set()
    print("  verified: example CAN be repaired at n=3 with irreducible symmetric c.")

    def tilde3(i, n):
        if (i, n) == (0, 2):
            return canonical(0, 2) ^ bd_zeta1
        if (i, n) == (0, 3):
            return canonical(0, 3) ^ sol
        return canonical(i, n)

    bad = []
    for n in range(NMAX + 1):
        for i in range(NMAX + 1):
            d = check_identity(tilde3, i, n)
            if d:
                bad.append((i, n, len(d)))
    print(f"  with this repair, failing identities: {bad if bad else 'none up to n=%d' % NMAX}")
    t03 = tilde3(0, 3)
    print(f"  repaired Delta_0[3]: irreducible={is_irreducible(t03)}, free={is_free(t03, 0, 3)}")
