"""Check the combinatorial claims in the simpliciality proof, and a proposed
reformulation via the transposition automorphism tau exchanging x_j, x_{j+1}."""
from itertools import combinations

def subsets(n):
    return [frozenset(S) for r in range(n + 2) for S in combinations(range(n + 1), r)]

def codeg(U, j, n):
    """P(sigma_j)(U) as a subset of {0..n-1}, or None if zero."""
    if j + 1 in U:
        keep = U - {j + 1}
    elif j in U:
        keep = U - {j}
    else:
        return None
    return frozenset(u if u <= j else u - 1 for u in keep)

def tau(U, j):
    """exchange j and j+1"""
    out = set(U)
    a, b = j in U, j + 1 in U
    out.discard(j); out.discard(j + 1)
    if a: out.add(j + 1)
    if b: out.add(j)
    return frozenset(out)

for n in range(1, 6):
    S = subsets(n)
    # (a) P(sigma_j) o tau = P(sigma_j)
    a_ok = all(codeg(tau(U, j), j, n) == codeg(U, j, n)
               for j in range(n) for U in S)
    # (b) nonzero-image irreducibles come in tau-pairs; exactly two per image
    b_ok = True
    for j in range(n):
        img = {}
        for V in S:
            for W in S:
                if V & W:
                    continue
                cv, cw = codeg(V, j, n), codeg(W, j, n)
                if cv is None or cw is None:
                    continue
                img.setdefault((cv, cw), []).append((V, W))
        for key, pre in img.items():
            if len(pre) != 2:
                b_ok = False
            V, W = pre[0]
            if (tau(V, j), tau(W, j)) != pre[1]:
                b_ok = False
            # each such element has one of j,j+1 in V and the other in W
            if not ({j, j + 1} & V and {j, j + 1} & W):
                b_ok = False
    print(f"n={n}: P(sigma_j).tau = P(sigma_j): {a_ok};  "
          f"exactly two preimages, swapped by tau: {b_ok}")
