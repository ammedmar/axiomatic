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


