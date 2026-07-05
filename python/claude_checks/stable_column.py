"""Stable single-column solutions: image of restriction from N=6 (and 7) down.
Extract a non-twist generator at i=0 and display it."""
from single_column import column_space, irr_deg
from cupi import one_plus_T, canonical

def restrict_rank(i, M, N, return_span=False):
    unknowns, idx, dim, basis = column_space(i, M)
    mask = 0
    for k, (n, V, W) in enumerate(unknowns):
        if n <= N:
            mask |= 1 << k
    restricted = [v & mask for v in basis]
    pivots = {}
    for r in restricted:
        while r:
            p = r.bit_length() - 1
            if p in pivots: r ^= pivots[p]
            else: pivots[p] = r; break
    if return_span:
        return len(pivots), pivots, unknowns, idx, mask
    return len(pivots)

for i in range(3):
    r64 = restrict_rank(i, 6, 4)
    r54 = restrict_rank(i, 5, 4)
    print(f"i={i}: dim at N=4: {column_space(i,4)[2]}, image from N=5: {r54}, image from N=6: {r64}")
