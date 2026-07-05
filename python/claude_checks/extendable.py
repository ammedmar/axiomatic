"""For each pair N < M: dimension of the image of the restriction map
{irreducible constructions truncated at M} -> {level <= N data}.
If this stabilizes at 2 (= span{canonical, T canonical}), then every
irreducible construction agrees with a linear combination of canonical and
T canonical in low degrees -- evidence that non-degenerate + irreducible
implies uniqueness (freeness removable under the stronger axiom)."""
from solution_space import solve, construction_to_vec, irr_unknowns
from cupi import canonical, T

MMAX = 6

for M in range(3, MMAX + 1):
    unknowns, idx, dim_null, null_basis = solve(M)
    print(f"level M={M}: dim(solutions)={dim_null}")
    for N in range(2, M):
        # mask of unknowns with n <= N
        mask = 0
        for k, (n, V, W) in enumerate(unknowns):
            if n <= N:
                mask |= 1 << k
        restricted = [v & mask for v in null_basis]
        # rank over F_2
        pivots = {}
        for r in restricted:
            while r:
                p = r.bit_length() - 1
                if p in pivots:
                    r ^= pivots[p]
                else:
                    pivots[p] = r
                    break
        rank = len(pivots)
        # is span == span{can, Tcan} restricted? check both reduce to 0
        can_vec = construction_to_vec(canonical, N, {u: k for k, u in idx.items() if False} or idx) & mask
        tcan_vec = construction_to_vec(lambda i, n: T(canonical(i, n)), N, idx) & mask
        def reduces_to_zero(v):
            while v:
                p = v.bit_length() - 1
                if p not in pivots:
                    return False
                v ^= pivots[p]
            return True
        in_span = reduces_to_zero(can_vec) and reduces_to_zero(tcan_vec)
        print(f"  restriction to N={N}: dim(image)={rank}, contains can & Tcan: {in_span}")
