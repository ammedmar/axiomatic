"""For X[n] = zeta_{c(u)} + zeta_{c(v)} (co-vertex pair, fresh at level n,
single column i=0): which pairs (u,v) extend 1, 2, 3 levels up?"""
from itertools import combinations
from cupi import partition_chain, cofaces_sum, one_plus_T
from single_column import irr_deg

def covertex_zeta(u, n):
    U = frozenset(range(n + 1)) - {u}
    return partition_chain(U)

def solve_affine(fixed, n_start, N, i0=0):
    """fixed: dict level->chain for levels <= n_start. Unknowns: levels n_start+1..N.
    Single-column system at index i0. Returns (consistent, particular dict, freedom dim)."""
    unknowns = []
    for n in range(n_start + 1, N + 1):
        for b in irr_deg(n, i0):
            unknowns.append((n, b[0], b[1]))
    idx = {u: k for k, u in enumerate(unknowns)}
    eqs, rhs = {}, {}
    def touch(key): return eqs.setdefault(key, 0)
    for k, (n, V, W) in enumerate(unknowns):
        eqs[('s', n, V, W)] = touch(('s', n, V, W)) ^ (1 << k)
        eqs[('s', n, W, V)] = touch(('s', n, W, V)) ^ (1 << k)
        if len(V) < n:
            for u in range(n + 1):
                if u not in V:
                    key = ('c', n, frozenset(V | {u}), W); eqs[key] = touch(key) ^ (1 << k)
        if len(W) < n:
            for u in range(n + 1):
                if u not in W:
                    key = ('c', n, V, frozenset(W | {u})); eqs[key] = touch(key) ^ (1 << k)
        if n + 1 <= N:
            for j in range(n + 2):
                V2 = frozenset({x if x < j else x + 1 for x in V} | {j})
                W2 = frozenset({x if x < j else x + 1 for x in W} | {j})
                key = ('c', n + 1, V2, W2); eqs[key] = touch(key) ^ (1 << k)
    # inhomogeneous: cofaces of fixed level n_start into level n_start+1
    ch = fixed.get(n_start, set())
    if ch and n_start + 1 <= N:
        for (V, W) in cofaces_sum(ch, n_start + 1):
            key = ('c', n_start + 1, V, W)
            touch(key)
            rhs[key] = rhs.get(key, 0) ^ 1
    rows = [(m, rhs.get(key, 0)) for key, m in eqs.items() if m or rhs.get(key, 0)]
    pivots = {}
    for m, b in rows:
        while m:
            p = m.bit_length() - 1
            if p in pivots:
                pm, pb = pivots[p]; m ^= pm; b ^= pb
            else:
                pivots[p] = (m, b); break
        else:
            if b:
                return False, None, None
    part = 0
    for p, (m, b) in sorted(pivots.items()):
        if b ^ (bin((m ^ (1 << p)) & part).count('1') % 2):
            part |= 1 << p
    comps = {}
    for k, (n, V, W) in enumerate(unknowns):
        if part >> k & 1:
            comps.setdefault(n, set()).add((V, W))
    return True, comps, len(unknowns) - len(pivots)

# which co-vertex pairs at level n extend how far?
for n in (3, 4, 5):
    print(f"level n={n}, pairs (u,v), extendable to N=n+1 / n+2:")
    row = []
    for u in range(n + 1):
        for v in range(u + 1, n + 1):
            X = covertex_zeta(u, n) ^ covertex_zeta(v, n)
            ok1, _, _ = solve_affine({n: X}, n, n + 1)
            ok2, _, _ = solve_affine({n: X}, n, n + 2)
            row.append(f"({u},{v}):{'+' if ok1 else '-'}{'+' if ok2 else '-'}")
    print("  " + "  ".join(row))
