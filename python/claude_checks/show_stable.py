"""Display stable single-column generators at i=0, levels up to 5, from M=6."""
from single_column import column_space
from cupi import one_plus_T, canonical, partition_chain

M, N = 6, 5
i0 = 0
unknowns, idx, dim, basis = column_space(i0, M)
mask = 0
for k, (n, V, W) in enumerate(unknowns):
    if n <= N:
        mask |= 1 << k

# twist generator restricted
tw = 0
for n in range(N + 1):
    for (V, W) in one_plus_T(canonical(i0, n)):
        tw |= 1 << idx[(n, V, W)]

# span of restricted stable solutions
restricted = sorted({v & mask for v in basis} - {0})
pivots = {}
order = []
for r0 in restricted:
    r = r0
    while r:
        p = r.bit_length() - 1
        if p in pivots: r ^= pivots[p]
        else: pivots[p] = r; order.append(r); break
print(f"stable dim at N={N} from M={M}: {len(pivots)}")

# reduce twist against span to confirm membership
r = tw
while r:
    p = r.bit_length() - 1
    if p not in pivots: break
    r ^= pivots[p]
print(f"twist in span: {r == 0}")

# build quotient-by-twist representatives: reduce each spanning vector by tw where helpful
def comps_of(vec):
    comps = {}
    for k, (n, V, W) in enumerate(unknowns):
        if vec >> k & 1:
            comps.setdefault(n, set()).add((V, W))
    return comps

def show(vec, name):
    comps = comps_of(vec)
    print(f"=== {name} ===")
    for n in sorted(comps):
        terms = sorted(comps[n], key=lambda x: (len(x[0]), sorted(x[0]), sorted(x[1])))
        s = " + ".join(f"{tuple(sorted(V))}|{tuple(sorted(W))}" for V, W in terms)
        print(f"  X[{n}] ({len(terms)} terms): {s}")

# canonical basis of the stable span, then pick vectors with small low-level support
span_vecs = order
# generate all 2^4 combos, pick the ones minimizing support at low levels, exclude 0 and tw
from itertools import combinations as icomb
best = []
allv = []
for rbits in range(1, 1 << len(span_vecs)):
    v = 0
    for b in range(len(span_vecs)):
        if rbits >> b & 1:
            v ^= span_vecs[b]
    allv.append(v)
allv = sorted(set(allv) - {0, tw}, key=lambda v: bin(v).count('1'))
for v in allv[:3]:
    show(v, f"stable non-twist generator (weight {bin(v).count('1')})")
