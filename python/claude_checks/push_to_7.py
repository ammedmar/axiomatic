"""Extend the fresh-at-3 pair and the all-co-vertices sum to N=7, single column i=0."""
from pair_extend import covertex_zeta, solve_affine

n = 3
pair = covertex_zeta(0, n) ^ covertex_zeta(1, n)
allc = covertex_zeta(0, n) ^ covertex_zeta(1, n) ^ covertex_zeta(2, n) ^ covertex_zeta(3, n)

for name, X in (("pair z_c(0)+z_c(1)", pair), ("all co-vertices A[3]", allc)):
    for N in (5, 6, 7):
        ok, comps, freedom = solve_affine({n: X}, n, N)
        print(f"{name}: extend 3->{N}: consistent={ok}" + (f", freedom={freedom}" if ok else ""))
    print()
