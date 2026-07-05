"""Verify the F-language reformulation of single-column (i=0) solutions.

X[n] = sum of F_n(u; V,W) * (V ot W) over partitions V|W of c(u) := {0..n}\{u}.

Claimed equivalence of the defining identities with:
(R) recurrence: for k in c(u), partitions (A,B) of c(u)\{k}:
    F_n(u; A+k, B) + F_n(u; A, B+k) = F_{n-1}(s_k(u); s_k(A), s_k(B))
(C) cycle: for every full partition (P,Q) of {0..n}, P,Q nonempty:
    sum_{u in P} F_n(u; P-u, Q) + sum_{u in Q} F_n(u; P, Q-u) = 0
(S) symmetry: F_n(u; V, W) = F_n(u; W, V).
"""
from itertools import combinations
from pair_extend import covertex_zeta, solve_affine

def sk(S, k):
    return frozenset(x if x < k else x - 1 for x in S)

n0 = 3
X3 = covertex_zeta(0, n0) ^ covertex_zeta(1, n0)
ok, comps, _ = solve_affine({n0: X3}, n0, 5)
comps[3] = X3
F = {}  # (n, u, V, W) -> 1 if present
for n, ch in comps.items():
    for (V, W) in ch:
        u = next(iter(frozenset(range(n + 1)) - V - W))
        F[(n, u, V, W)] = 1

def getF(n, u, V, W):
    return F.get((n, u, V, W), 0)

# check (S), (R), (C) for n = 4, 5
allok = True
for n in (4, 5):
    for u in range(n + 1):
        cu = [x for x in range(n + 1) if x != u]
        for r in range(len(cu) + 1):
            for V in combinations(cu, r):
                V = frozenset(V); W = frozenset(cu) - V
                if getF(n, u, V, W) != getF(n, u, W, V):
                    allok = False; print(f"(S) fails n={n} u={u}")
        for k in cu:
            rest = [x for x in cu if x != k]
            for r in range(len(rest) + 1):
                for A in combinations(rest, r):
                    A = frozenset(A); B = frozenset(rest) - A
                    lhs = getF(n, u, A | {k}, B) ^ getF(n, u, A, B | {k})
                    u2 = u if u < k else u - 1
                    rhs = getF(n - 1, u2, sk(A, k), sk(B, k))
                    if lhs != rhs:
                        allok = False; print(f"(R) fails n={n} u={u} k={k} A={sorted(A)}")
    # cycle condition
    full = list(range(n + 1))
    for r in range(1, n + 1):
        for P in combinations(full, r):
            P = frozenset(P); Q = frozenset(full) - P
            s = 0
            for u in P:
                s ^= getF(n, u, P - {u}, Q)
            for u in Q:
                s ^= getF(n, u, P, Q - {u})
            if s:
                allok = False; print(f"(C) fails n={n} P={sorted(P)}")
print("F-language conditions (S),(R),(C) verified on computed solution:", allok)
