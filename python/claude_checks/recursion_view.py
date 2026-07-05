"""Display Y[4], Y[5] above X[3] = zeta_{c(0)} + zeta_{c(1)}, canonicalized."""
from pair_extend import covertex_zeta, solve_affine
from cupi import partition_chain, one_plus_T, T

n = 3
X = covertex_zeta(0, n) ^ covertex_zeta(1, n)
ok, comps, freedom = solve_affine({n: X}, n, n + 2)
print(f"extend {n}->{n+2}: consistent={ok}, freedom dim={freedom}")

def decompose(ch):
    """Split chain into full partition chains + symmetric rest + asymmetric rest."""
    rest = set(ch); pcs = []
    changed = True
    while changed:
        changed = False
        for (V, W) in list(rest):
            zU = partition_chain(V | W)
            if zU <= rest:
                rest -= zU; pcs.append(tuple(sorted(V | W))); changed = True
    sym = {b for b in rest if (b[1], b[0]) in rest}
    asym = rest - sym
    return pcs, sym, asym

for lvl in sorted(comps):
    ch = comps[lvl]
    pcs, sym, asym = decompose(ch)
    print(f"\nY[{lvl}]: {len(ch)} terms, symmetric: {one_plus_T(ch) == set()}")
    if pcs:
        print(f"  full partition chains: {' + '.join('z' + str(p) for p in sorted(pcs))}")
    if sym:
        bysupp = {}
        for V, W in sym:
            bysupp.setdefault(tuple(sorted(V | W)), []).append((tuple(sorted(V)), tuple(sorted(W))))
        for supp in sorted(bysupp):
            terms = sorted(bysupp[supp])
            print(f"  sym rest on U={supp} ({len(terms)}/{2**len(supp)}): {' + '.join(f'{a}|{b}' for a, b in terms)}")
    if asym:
        print(f"  ASYMMETRIC rest: {len(asym)} terms")
