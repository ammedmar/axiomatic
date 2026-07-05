"""Solve the single-column problem entirely in constants space.

Unknowns: c_m(w), 0 <= w <= m <= N.
Rule reconstruction: F_n(u)(V) = sum_{A subset V} c_{n-|A|}(u - |A cap [0,u)|),
type-only: f(n,u,p,q) = sum_{j<=p,l<=q} C(p,j)C(q,l) c_{n-j-l}(u-j).

Conditions:
(S) for each n <= N, u <= n: sum_{(j,l) != (0,0)} C(u,j)C(n-u,l) c_{n-j-l}(u-j) = 0
(C) for each n <= N: Phi_n([0,n),{n}) = sum_{u<n} f(n,u,u,n-1-u) + f(n,n,n,0) = 0
"""
from math import comb

N = 12
idx = {}
for m in range(N + 1):
    for w in range(m + 1):
        idx[(m, w)] = len(idx)
NV = len(idx)

def f_row(n, u, p, q):
    """Linear functional (bitmask over constants) for f(n,u,p,q)."""
    row = 0
    for j in range(p + 1):
        if comb(p, j) % 2 == 0: continue
        for l in range(q + 1):
            if comb(q, l) % 2 == 0: continue
            m, w = n - j - l, u - j
            if 0 <= w <= m:
                row ^= 1 << idx[(m, w)]
    return row

rows = []
for n in range(N + 1):
    # (S)
    for u in range(n + 1):
        row = f_row(n, u, u, n - u) ^ (1 << idx[(n, u)])  # full transform minus (0,0) term
        if row:
            rows.append(row)
    # (C)
    if n >= 1:
        row = 0
        for u in range(n):
            row ^= f_row(n, u, u, n - 1 - u)
        row ^= f_row(n, n, n, 0)
        if row:
            rows.append(row)

pivots = {}
for r in rows:
    while r:
        p = r.bit_length() - 1
        if p in pivots: r ^= pivots[p]
        else: pivots[p] = r; break
dim = NV - len(pivots)
print(f"N={N}: constants unknowns={NV}, dim(solution space)={dim}")

free_vars = [k for k in range(NV) if k not in pivots]
basis = []
for fv in free_vars:
    vec = 1 << fv
    for p, r in sorted(pivots.items()):
        if bin((r ^ (1 << p)) & vec).count('1') % 2:
            vec |= 1 << p
    basis.append(vec)
for vec in basis:
    assert all(bin(r & vec).count('1') % 2 == 0 for r in rows)

# twist check: c_m(w) = [w=0]+[w=m] (m>=1), c_0(0)=0
tw = 0
for m in range(1, N + 1):
    for w in (0, m):
        tw ^= 1 << idx[(m, w)]
def in_span(v):
    for p, r in sorted(pivots.items(), reverse=True):
        pass
    x = v
    while x:
        p = x.bit_length() - 1
        # reduce against basis span: build span pivots
        return None
# simpler: check twist satisfies all rows
print("twist satisfies all conditions:", all(bin(r & tw).count('1') % 2 == 0 for r in rows))

# find solution with c_3 = (1,1,0,0), c_m = 0 for m < 3: add those as constraints
extra = []
for m in range(3):
    for w in range(m + 1):
        extra.append((1 << idx[(m, w)], 0))
extra.append((1 << idx[(3, 0)], 1))
extra.append((1 << idx[(3, 1)], 1))
extra.append((1 << idx[(3, 2)], 0))
extra.append((1 << idx[(3, 3)], 0))
arows = [(r, 0) for r in rows] + extra
piv = {}
ok = True
for m, b in arows:
    while m:
        p = m.bit_length() - 1
        if p in piv:
            pm, pb = piv[p]; m ^= pm; b ^= pb
        else:
            piv[p] = (m, b); break
    else:
        if b: ok = False
print(f"pair family c_3=(1,1,0,0) extends in constants space to N={N}:", ok)
if ok:
    part = 0
    for p, (m, b) in sorted(piv.items()):
        if b ^ (bin((m ^ (1 << p)) & part).count('1') % 2):
            part |= 1 << p
    # display triangle
    print("one valid triangle (rows m=3..N, entries c_m(0..m)):")
    for m in range(3, N + 1):
        vals = [(part >> idx[(m, w)]) & 1 for w in range(m + 1)]
        print(f"  m={m:2d}: {''.join(map(str, vals))}")
