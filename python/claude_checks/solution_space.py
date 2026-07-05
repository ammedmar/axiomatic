"""
1. Check (1+T)canonical is a valid semi-simplicial cup-i construction:
   irreducible, non-zero, NOT free, but degenerate (x cup_0 x = 0 on vertices).

2. Compute the F_2-dimension of the space of ALL irreducible truncated
   (n <= N) semi-simplicial cup-i constructions. The identities are linear
   and irreducibility is a linear (support) condition, so solutions form a
   vector space. If its dimension is 2 = dim span{canonical, T canonical},
   then non-degenerate + irreducible uniquely determine the construction up
   to iso at truncation level N (freeness not needed once non-degeneracy is
   imposed). If dim > 2, inspect extra solutions.
"""

from itertools import combinations
from cupi import (bd, bd_subset, T, one_plus_T, coface_subset,
                           cofaces_sum, canonical, check_identity,
                           is_irreducible, is_free)

NMAX = 5

print("=== A. (1+T)canonical as a construction ===")


def opt(i, n):
    return one_plus_T(canonical(i, n))


bad = []
for n in range(NMAX + 1):
    for i in range(NMAX + 2):
        d = check_identity(opt, i, n)
        if d:
            bad.append((i, n))
print(f"  identities up to n={NMAX}: {'all pass' if not bad else 'FAIL at %s' % bad}")
nonzero = any(opt(i, n) for n in range(NMAX + 1) for i in range(n + 1))
print(f"  non-zero: {nonzero}")
print(f"  irreducible everywhere: "
      f"{all(is_irreducible(opt(i, n)) for n in range(NMAX + 1) for i in range(n + 1))}")
notfree = [(i, n) for n in range(NMAX + 1) for i in range(n + 1)
           if not is_free(opt(i, n), i, n)]
print(f"  freeness fails at (i,n): {notfree[:5]}{' ...' if len(notfree) > 5 else ''}")
print(f"  Delta_0[0] = {opt(0, 0)}  (empty => degenerate: x cup_0 x = 0 on vertices)")
print(f"  simplicial check via sigma_j is skipped (not needed for the counterexample)")

print()
print("=== B. dimension of truncated irreducible construction space ===")


def irr_unknowns(N):
    """All (n, V, W): V, W disjoint subsets of {0..n}, |V|+|W| <= n.
    (i = n - |V| - |W| is determined by the triple.)"""
    out = []
    for n in range(N + 1):
        elts = list(range(n + 1))
        for kv in range(n + 2):
            for V in combinations(elts, kv):
                rest = [u for u in elts if u not in V]
                for kw in range(n - kv + 1):
                    for W in combinations(rest, kw):
                        out.append((n, frozenset(V), frozenset(W)))
    return out


def solve(N, extra_rows=None, inhomog=None):
    """Build the linear system for truncated constructions supported on
    irreducible basis elements; return (nullspace basis as lists of unknown
    indices) via F_2 Gaussian elimination with bitmask rows over unknowns."""
    unknowns = irr_unknowns(N)
    idx = {u: k for k, u in enumerate(unknowns)}
    nun = len(unknowns)

    # equations: dict (n, i, targetV, targetW) -> bitmask of unknowns
    eqs = {}

    def touch(key):
        return eqs.setdefault(key, 0)

    for k, (n, V, W) in enumerate(unknowns):
        i = n - len(V) - len(W)
        # contribution to identity (i, n) via bd Delta_i[n]
        for V2 in bd_subset(V, n):
            key = (n, i, V2, W)
            eqs[key] = touch(key) ^ (1 << k)
        for W2 in bd_subset(W, n):
            key = (n, i, V, W2)
            eqs[key] = touch(key) ^ (1 << k)
        # contribution to identity (i, n+1) via sum_j P(delta_j)^{ot2} Delta_i[n]
        if n + 1 <= N:
            for j in range(n + 2):
                key = (n + 1, i, coface_subset(V, j), coface_subset(W, j))
                eqs[key] = touch(key) ^ (1 << k)
        # contribution to identity (i+1, n) via (1+T) Delta_i[n]
        eqs[(n, i + 1, V, W)] = touch((n, i + 1, V, W)) ^ (1 << k)
        eqs[(n, i + 1, W, V)] = touch((n, i + 1, W, V)) ^ (1 << k)

    rows = [r for r in eqs.values() if r]

    # Gaussian elimination, track pivots
    pivot_of = {}  # bit position -> reduced row
    for r in rows:
        while r:
            p = r.bit_length() - 1
            if p in pivot_of:
                r ^= pivot_of[p]
            else:
                pivot_of[p] = r
                break
    rank = len(pivot_of)
    dim_null = nun - rank

    # nullspace basis: free variables = non-pivot columns
    free_vars = [k for k in range(nun) if k not in pivot_of]
    null_basis = []
    # back-substitution: each pivot row has its pivot as HIGHEST bit, so the
    # pivot value depends only on lower bits -> process pivots in ASCENDING order
    sorted_pivots = sorted(pivot_of.items())
    for fv in free_vars:
        vec = 1 << fv
        for p, r in sorted_pivots:
            # value of pivot var p: parity of (r without pivot bit) . vec
            if bin((r ^ (1 << p)) & vec).count('1') % 2:
                vec |= 1 << p
        null_basis.append(vec)
    # hard verification: every null basis vector satisfies every equation
    all_rows = [r for r in eqs.values() if r]
    for vec in null_basis:
        assert all(bin(r & vec).count('1') % 2 == 0 for r in all_rows), \
            "nullspace vector fails an equation"
    return unknowns, idx, dim_null, null_basis


def construction_to_vec(D, N, idx):
    vec = 0
    for n in range(N + 1):
        for i in range(n + 1):
            for (V, W) in D(i, n):
                vec |= 1 << idx[(n, V, W)]
    return vec


for N in range(2, NMAX + 1):
    unknowns, idx, dim_null, null_basis = solve(N)
    can_vec = construction_to_vec(canonical, N, idx)
    tcan_vec = construction_to_vec(lambda i, n: T(canonical(i, n)), N, idx)
    # check canonical and T canonical are in the nullspace span
    # (they satisfy the equations by construction; verify via reduction)
    print(f"  N={N}: unknowns={len(unknowns)}, dim(solution space)={dim_null}")
    if dim_null <= 6:
        # describe each basis solution: where does it deviate from 0, degenerate?
        for b, vec in enumerate(null_basis):
            supp = [unknowns[k] for k in range(len(unknowns)) if vec >> k & 1]
            deg00 = (0, frozenset(), frozenset())
            nondeg = idx[deg00] if deg00 in idx else None
            has00 = bool(vec >> idx[deg00] & 1)
            eq_can = vec == can_vec
            eq_tcan = vec == tcan_vec
            eq_sum = vec == can_vec ^ tcan_vec
            tag = ("= canonical" if eq_can else
                   "= T canonical" if eq_tcan else
                   "= (1+T) canonical" if eq_sum else
                   f"OTHER ({len(supp)} terms)")
            print(f"    basis vector {b}: nondegenerate={has00}, {tag}")
