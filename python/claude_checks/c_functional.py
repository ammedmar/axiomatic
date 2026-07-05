"""Extract the (C)_n functional coefficients on constants c_m(w); guess pattern."""
from math import comb

N = 12
idx = {}
pos = []
for m in range(N + 1):
    for w in range(m + 1):
        idx[(m, w)] = len(idx); pos.append((m, w))

def f_row(n, u, p, q):
    row = 0
    for j in range(p + 1):
        if comb(p, j) % 2 == 0: continue
        for l in range(q + 1):
            if comb(q, l) % 2 == 0: continue
            m, w = n - j - l, u - j
            if 0 <= w <= m:
                row ^= 1 << idx[(m, w)]
    return row

for n in range(2, 11):
    row = 0
    for u in range(n):
        row ^= f_row(n, u, u, n - 1 - u)
    row ^= f_row(n, n, n, 0)
    coeffs = [(m, w) for k, (m, w) in enumerate(pos) if row >> k & 1]
    print(f"(C)_{n}: " + " ".join(f"c_{m}({w})" for m, w in sorted(coeffs)))
