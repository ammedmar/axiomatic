"""Extend a truncated irreducible construction one level up by solving the
affine system at level n. Canonicalize modulo the level-n kernel."""
from itertools import combinations
from cupi import bd, T, one_plus_T, cofaces_sum, partition_chain, canonical, check_identity
from kernel_check import irr_level, level_only_solutions

def solve_level(D, n):
    """D: dict (i,m) -> chain for m < n. Returns (particular, nullbasis, unknowns)."""
    unknowns = irr_level(n)
    idx = {u: k for k, u in enumerate(unknowns)}
    nun = len(unknowns)
    eqs = {}   # key -> bitmask over unknowns
    rhs = {}   # key -> F2 constant
    def touch(key): return eqs.setdefault(key, 0)
    for k, (V, W) in enumerate(unknowns):
        i = n - len(V) - len(W)
        if len(V) < n:
            for u in range(n + 1):
                if u not in V:
                    key = (i, frozenset(V | {u}), W); eqs[key] = touch(key) ^ (1 << k)
        if len(W) < n:
            for u in range(n + 1):
                if u not in W:
                    key = (i, V, frozenset(W | {u})); eqs[key] = touch(key) ^ (1 << k)
        eqs[(i+1, V, W)] = touch((i+1, V, W)) ^ (1 << k)
        eqs[(i+1, W, V)] = touch((i+1, W, V)) ^ (1 << k)
    # constants: coface pushforwards of D(i, n-1)
    for i in range(n + 2):
        ch = D.get((i, n - 1), set())
        if ch:
            push = cofaces_sum(ch, n)
            for (V, W) in push:
                key = (i, V, W)
                touch(key)
                rhs[key] = rhs.get(key, 0) ^ 1
    # Gaussian elimination on [A | b]
    rows = [(mask, rhs.get(key, 0)) for key, mask in eqs.items() if mask or rhs.get(key, 0)]
    pivots = {}
    for mask, b in rows:
        while mask:
            p = mask.bit_length() - 1
            if p in pivots:
                pm, pb = pivots[p]
                mask ^= pm; b ^= pb
            else:
                pivots[p] = (mask, b); break
        else:
            if b:
                return None, None, unknowns  # inconsistent
    # particular solution: free vars = 0
    part = 0
    for p, (mask, b) in sorted(pivots.items()):
        val = b ^ (bin((mask ^ (1 << p)) & part).count('1') % 2)
        if val:
            part |= 1 << p
    # verify
    for key, mask in eqs.items():
        assert bin(mask & part).count('1') % 2 == rhs.get(key, 0), "particular fails"
    _, _, dim, nullb = level_only_solutions(n)
    return part, nullb, unknowns

def vec_to_comps(vec, unknowns, n):
    comps = {}
    for k, (V, W) in enumerate(unknowns):
        if vec >> k & 1:
            i = n - len(V) - len(W)
            comps.setdefault(i, set()).add((V, W))
    return comps

def minimize(part, nullb):
    """Greedy: reduce particular support using nullspace vectors."""
    best = part
    improved = True
    while improved:
        improved = False
        for v in nullb:
            cand = best ^ v
            if bin(cand).count('1') < bin(best).count('1'):
                best = cand; improved = True
    return best

def show(comps, n, label):
    print(f"--- {label} (level {n}) ---")
    for i in sorted(comps):
        ch = comps[i]
        # try to express: split off complete partition chains
        rest = set(ch); pcs = []
        changed = True
        while changed:
            changed = False
            for (V, W) in list(rest):
                U = V | W
                zU = partition_chain(U)
                zU = {b for b in zU}  # all valid here
                if zU <= rest:
                    rest -= zU; pcs.append(tuple(sorted(U))); changed = True
        pcs_str = " + ".join(f"z{p}" for p in sorted(pcs)) if pcs else ""
        rest_str = " + ".join(f"{tuple(sorted(V))}|{tuple(sorted(W))}"
                              for V, W in sorted(rest, key=lambda x: (sorted(x[0]), sorted(x[1]))))
        joiner = " + " if pcs_str and rest_str else ""
        print(f"  Y_{i} [{len(ch)} terms] = {pcs_str}{joiner}{rest_str}")

# ---- the generator G at level 2: X_0 = zeta_{01} + zeta_{02} ----
f = frozenset
G02 = partition_chain({0, 1}) ^ partition_chain({0, 2})
D = {(0, 2): G02}

part, nullb, unknowns = solve_level(D, 3)
if part is None:
    print("LEVEL 3: INCONSISTENT — no extension exists!")
else:
    m = minimize(part, nullb)
    comps3 = vec_to_comps(m, unknowns, 3)
    show(comps3, 3, "minimal extension of G = bd zeta_{0} at (0,2)")
    # store and go one more level
    for i, ch in comps3.items():
        D[(i, 3)] = ch
    part4, nullb4, unknowns4 = solve_level(D, 4)
    if part4 is None:
        print("LEVEL 4: INCONSISTENT with this level-3 choice")
    else:
        m4 = minimize(part4, nullb4)
        comps4 = vec_to_comps(m4, unknowns4, 4)
        show(comps4, 4, "minimal extension continued")
