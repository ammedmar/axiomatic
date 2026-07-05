"""Impose the simplicial constraint P(sigma_j)^{ot2} Delta_i[n] = 0 as extra
linear equations and recompute solution space dimensions."""
from itertools import combinations
from solution_space import solve, construction_to_vec, irr_unknowns
from cupi import canonical, T, one_plus_T, bd_subset, coface_subset, is_free

def codegeneracy_subset(U, j):
    """P(sigma_j): returns subset of {0..n-1} or None (=0)."""
    if j + 1 in U:
        keep = U - {j + 1}
    elif j in U:
        keep = U - {j}
    else:
        return None
    return frozenset(u if u <= j else u - 1 for u in keep)

# sanity: canonical must be simplicial
ok = True
for n in range(1, 7):
    for i in range(n + 1):
        for j in range(n):
            out = set()
            for V, W in canonical(i, n):
                V2 = codegeneracy_subset(V, j)
                W2 = codegeneracy_subset(W, j)
                if V2 is not None and W2 is not None:
                    out ^= {(V2, W2)}
            if out:
                ok = False
                print(f"canonical fails simpliciality at (i={i},n={n},j={j})")
print(f"sanity: canonical simplicial up to n=6: {ok}")

def solve_simplicial(M):
    unknowns = irr_unknowns(M)
    idx = {u: k for k, u in enumerate(unknowns)}
    eqs = {}
    def touch(key):
        return eqs.setdefault(key, 0)
    for k, (n, V, W) in enumerate(unknowns):
        i = n - len(V) - len(W)
        for V2 in bd_subset(V, n):
            key = ('id', n, i, V2, W); eqs[key] = touch(key) ^ (1 << k)
        for W2 in bd_subset(W, n):
            key = ('id', n, i, V, W2); eqs[key] = touch(key) ^ (1 << k)
        if n + 1 <= M:
            for j in range(n + 2):
                key = ('id', n + 1, i, coface_subset(V, j), coface_subset(W, j))
                eqs[key] = touch(key) ^ (1 << k)
        eqs[('id', n, i + 1, V, W)] = touch(('id', n, i + 1, V, W)) ^ (1 << k)
        eqs[('id', n, i + 1, W, V)] = touch(('id', n, i + 1, W, V)) ^ (1 << k)
        # simplicial equations
        for j in range(n):
            V2 = codegeneracy_subset(V, j)
            W2 = codegeneracy_subset(W, j)
            if V2 is not None and W2 is not None:
                key = ('deg', n, i, j, V2, W2)
                eqs[key] = touch(key) ^ (1 << k)
    rows = [r for r in eqs.values() if r]
    pivots = {}
    for r in rows:
        while r:
            p = r.bit_length() - 1
            if p in pivots:
                r ^= pivots[p]
            else:
                pivots[p] = r
                break
    nun = len(unknowns)
    dim_null = nun - len(pivots)
    free_vars = [k for k in range(nun) if k not in pivots]
    null_basis = []
    for fv in free_vars:
        vec = 1 << fv
        for p, r in sorted(pivots.items()):
            if bin((r ^ (1 << p)) & vec).count('1') % 2:
                vec |= 1 << p
        null_basis.append(vec)
    for vec in null_basis:
        assert all(bin(r & vec).count('1') % 2 == 0 for r in rows)
    return unknowns, idx, dim_null, null_basis

for M in range(2, 7):
    unknowns, idx, dim_null, null_basis = solve_simplicial(M)
    # twist family dim: canonical + w_i's that are nonzero at this level
    can_vec = construction_to_vec(canonical, M, idx)
    tw = [can_vec]
    for i0 in range(M + 1):
        vec = 0
        for n in range(M + 1):
            for (V, W) in one_plus_T(canonical(i0, n)):
                vec |= 1 << idx[(n, V, W)]
        if vec:
            tw.append(vec)
    # rank of twist span
    pv = {}
    for r in tw:
        while r:
            p = r.bit_length() - 1
            if p in pv:
                r ^= pv[p]
            else:
                pv[p] = r
                break
    print(f"M={M}: dim(simplicial solutions)={dim_null}, twist-family span={len(pv)}")
