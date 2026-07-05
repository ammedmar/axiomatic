"""Clean re-check: is there an irreducible c (symmetric) with bd c = defect D?"""
from itertools import combinations
from cupi import bd, T, one_plus_T, canonical, check_identity, partition_chain

zeta1 = partition_chain({1})
bd_zeta1 = bd(zeta1, 2)

def tilde(i, n):
    if (i, n) == (0, 2):
        return canonical(0, 2) ^ bd_zeta1
    return canonical(i, n)

D = check_identity(tilde, 0, 3)

# unknowns: symmetric orbits of irreducible degree-3 pairs over simplex^3
irr = []
for r in range(4):
    for V in combinations(range(4), r):
        rest = [u for u in range(4) if u not in V]
        for W in combinations(rest, 3 - r):
            irr.append((frozenset(V), frozenset(W)))
orbits, seen = [], set()
for b in irr:
    if b in seen: continue
    tb = (b[1], b[0]); seen |= {b, tb}
    orbits.append({b, tb})

# target space basis indexing
targets = {}
def tidx(x):
    return targets.setdefault(x, len(targets))

cols = []
for o in orbits:
    col = bd(set(o), 3)
    mask = 0
    for x in col:
        mask |= 1 << tidx(x)
    cols.append(mask)
Dmask = 0
for x in D:
    Dmask |= 1 << tidx(x)

# solve: is Dmask in span(cols)? proper elimination
pivots = {}
for r in cols:
    while r:
        p = r.bit_length() - 1
        if p in pivots: r ^= pivots[p]
        else: pivots[p] = r; break
v = Dmask
while v:
    p = v.bit_length() - 1
    if p not in pivots:
        print("CONFIRMED: no symmetric irreducible c with bd c = D (must also change Delta_1[3]).")
        break
    v ^= pivots[p]
else:
    print("solvable after all")
