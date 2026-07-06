"""Verification of sec/others.tex (final form).

1. Steenrod's construction (statement displays, interval k -> factor k mod 2)
   equals the canonical construction for all i <= n <= 8.
2. Real's construction with the corrected case split (BY PARITY OF i, not n):
   F1 = union of even-indexed gaps g_0, g_2, ...  (so g_i included iff i even),
   F2 = g_{-1} + union of odd-indexed gaps        (so g_i included iff i odd),
   equals the canonical construction for all i <= n <= 8.
3. Tie-break claims of both proofs.
"""
import sys
sys.path.insert(0, '/Users/laptop/Documents/research/active/axiomatic/python/claude_checks')
from itertools import combinations
from cupi import canonical, T

def steenrod(i, n):
    if i > n:
        return set()
    out = set()
    full = frozenset(range(n + 1))
    for ps in combinations(range(n + 1), i + 1):
        p = [0] + list(ps) + [n]
        f1, f2 = set(), set()
        for k in range(i + 2):
            (f1 if k % 2 == 0 else f2).update(range(p[k], p[k + 1] + 1))
        out ^= {(full - f1, full - f2)}
    return out

def real(i, n):
    if i > n:
        return set()
    out = set()
    for js in combinations(range(n + 1), i + 1):
        j = list(js)
        gaps = {-1: set(range(0, j[0]))}
        for k in range(i):
            gaps[k] = set(range(j[k] + 1, j[k + 1]))
        gaps[i] = set(range(j[i] + 1, n + 1))
        f1, f2 = set(), set(gaps[-1])
        for k in range(i + 1):
            (f1 if k % 2 == 0 else f2).update(gaps[k])
        out ^= {(frozenset(f1), frozenset(f2))}
    return out

okS = all(steenrod(i, n) == canonical(i, n) for n in range(9) for i in range(n + 1))
okR = all(real(i, n) == canonical(i, n) for n in range(9) for i in range(n + 1))
print(f"Steenrod == canonical for all i <= n <= 8: {okS}")
print(f"Real (i-parity split) == canonical for all i <= n <= 8: {okR}")

# tie-breaks: for every i <= 7, over simplex^{i+1}
f = frozenset
okT1 = okT2 = True
for i in range(8):
    n = i + 1
    tcan = T(canonical(i, n))
    # Steenrod proof: {0} ot empty is in T canonical but not in Steenrod's
    if (f({0}), f()) not in tcan or (f({0}), f()) in steenrod(i, n):
        okT1 = False
    # Real proof: empty ot {0} is in Real's but not in T canonical
    if (f(), f({0})) not in real(i, n) or (f(), f({0})) in tcan:
        okT2 = False
print(f"Steenrod tie-break ({{0}} ot empty in Tcan, not in Steenrod): {okT1}")
print(f"Real tie-break (empty ot {{0}} in Real, not in Tcan): {okT2}")
