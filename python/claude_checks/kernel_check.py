"""Verify: solutions supported only at level n = augmented partition-chain
cycles, dimension 2^n - 1."""
from itertools import combinations
from cupi import bd, T, one_plus_T, canonical, partition_chain, cofaces_sum

def irr_level(n):
    out = []
    for kv in range(n + 2):
        for V in combinations(range(n + 1), kv):
            rest = [u for u in range(n + 1) if u not in V]
            for kw in range(n + 1 - kv):
                for W in combinations(rest, kw):
                    out.append((frozenset(V), frozenset(W)))
    return out

def level_only_solutions(n):
    """Solve: bd X_i = (1+T) X_{i-1} for all i, X irreducible at level n only."""
    unknowns = irr_level(n)
    idx = {u: k for k, u in enumerate(unknowns)}
    eqs = {}
    def touch(key): return eqs.setdefault(key, 0)
    for k, (V, W) in enumerate(unknowns):
        i = n - len(V) - len(W)
        # bd X_i contributions to identity (i, n)
        for V2 in [frozenset(V | {u}) for u in range(n+1) if u not in V] if len(V) < n else []:
            key = (i, V2, W); eqs[key] = touch(key) ^ (1 << k)
        for W2 in [frozenset(W | {u}) for u in range(n+1) if u not in W] if len(W) < n else []:
            key = (i, V, W2); eqs[key] = touch(key) ^ (1 << k)
        # (1+T) X_i contributions to identity (i+1, n)
        eqs[(i+1, V, W)] = touch((i+1, V, W)) ^ (1 << k)
        eqs[(i+1, W, V)] = touch((i+1, W, V)) ^ (1 << k)
    rows = [r for r in eqs.values() if r]
    pivots = {}
    for r in rows:
        while r:
            p = r.bit_length() - 1
            if p in pivots: r ^= pivots[p]
            else: pivots[p] = r; break
    nun = len(unknowns)
    dim = nun - len(pivots)
    # nullspace basis
    free_vars = [k for k in range(nun) if k not in pivots]
    basis = []
    for fv in free_vars:
        vec = 1 << fv
        for p, r in sorted(pivots.items()):
            if bin((r ^ (1 << p)) & vec).count('1') % 2:
                vec |= 1 << p
        basis.append(vec)
    for vec in basis:
        assert all(bin(r & vec).count('1') % 2 == 0 for r in rows)
    return unknowns, idx, dim, basis

def is_partition_chain_sum(chain):
    """Check chain is a union of complete partition chains."""
    rest = set(chain)
    while rest:
        V, W = next(iter(rest))
        U = V | W
        zU = partition_chain(U)
        # truncate: keep only valid basis elements (parts of size <= n implicit: all are)
        if not zU <= rest:
            return False
        rest -= zU
    return True

for n in range(1, 6):
    unknowns, idx, dim, basis = level_only_solutions(n)
    all_pc = True
    for vec in basis:
        comps = {}
        for k, (V, W) in enumerate(unknowns):
            if vec >> k & 1:
                i = n - len(V) - len(W)
                comps.setdefault(i, set()).add((V, W))
        for i, ch in comps.items():
            if not is_partition_chain_sum(ch):
                all_pc = False
    print(f"n={n}: dim(level-only solutions)={dim}  (2^n-1={2**n-1})  all partition-chain sums: {all_pc}")
