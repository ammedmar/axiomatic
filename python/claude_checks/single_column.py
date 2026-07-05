"""Space of constructions supported at a single index i:
natural families X[n] in P(simplex^n)^{ot2}_{n+i}, irreducible support,
(1+T)X[n] = 0, and bd X[n] + sum_j P(delta_j)^{ot2} X[n-1] = 0.
Contains (1+T)canonical_i.  Question: anything else?"""
from itertools import combinations
from cupi import bd, T, one_plus_T, cofaces_sum, canonical

def irr_deg(n, i):
    """Irreducible basis elements of degree n+i over simplex^n: |V|+|W| = n-i."""
    out = []
    size = n - i
    if size < 0:
        return out
    for kv in range(size + 1):
        for V in combinations(range(n + 1), kv):
            rest = [u for u in range(n + 1) if u not in V]
            for W in combinations(rest, size - kv):
                out.append((frozenset(V), frozenset(W)))
    return out

def column_space(i, N):
    """Solve the single-column system for index i, levels n <= N."""
    unknowns = []
    for n in range(N + 1):
        for b in irr_deg(n, i):
            unknowns.append((n, b[0], b[1]))
    idx = {u: k for k, u in enumerate(unknowns)}
    eqs = {}
    def touch(key): return eqs.setdefault(key, 0)
    for k, (n, V, W) in enumerate(unknowns):
        # symmetry: (1+T)X[n] = 0
        eqs[('s', n, V, W)] = touch(('s', n, V, W)) ^ (1 << k)
        eqs[('s', n, W, V)] = touch(('s', n, W, V)) ^ (1 << k)
        # bd X[n] contributions to chain-map identity at level n
        if len(V) < n:
            for u in range(n + 1):
                if u not in V:
                    key = ('c', n, frozenset(V | {u}), W); eqs[key] = touch(key) ^ (1 << k)
        if len(W) < n:
            for u in range(n + 1):
                if u not in W:
                    key = ('c', n, V, frozenset(W | {u})); eqs[key] = touch(key) ^ (1 << k)
        # coface contributions to identity at level n+1
        if n + 1 <= N:
            for j in range(n + 2):
                V2 = frozenset({u if u < j else u + 1 for u in V} | {j})
                W2 = frozenset({u if u < j else u + 1 for u in W} | {j})
                key = ('c', n + 1, V2, W2); eqs[key] = touch(key) ^ (1 << k)
    rows = [r for r in eqs.values() if r]
    pivots = {}
    for r in rows:
        while r:
            p = r.bit_length() - 1
            if p in pivots: r ^= pivots[p]
            else: pivots[p] = r; break
    dim = len(unknowns) - len(pivots)
    # nullspace basis
    free_vars = [k for k in range(len(unknowns)) if k not in pivots]
    basis = []
    for fv in free_vars:
        vec = 1 << fv
        for p, r in sorted(pivots.items()):
            if bin((r ^ (1 << p)) & vec).count('1') % 2:
                vec |= 1 << p
        basis.append(vec)
    return unknowns, idx, dim, basis

for i in range(4):
    for N in (4, 5, 6):
        unknowns, idx, dim, basis = column_space(i, N)
        print(f"index i={i}, N={N}: dim(single-column solutions) = {dim}")
    print()
