"""Find an explicit non-degenerate irreducible construction (truncated at M)
outside the twist family {T^{eps_i} canonical}, and check its freeness."""
from solution_space import solve, construction_to_vec
from cupi import canonical, T, one_plus_T, is_free, is_irreducible

M = 6
unknowns, idx, dim_null, null_basis = solve(M)

can_vec = construction_to_vec(canonical, M, idx)

# twist family span: canonical and w_i = (1+T)canonical_i (single i component)
def w_vec(i0):
    vec = 0
    for n in range(M + 1):
        for (V, W) in one_plus_T(canonical(i0, n)):
            vec |= 1 << idx[(n, V, W)]
    return vec

trivial = [can_vec] + [w_vec(i) for i in range(M + 1)]
trivial = [v for v in trivial if v]

def make_reducer(vectors):
    pivots = {}
    for r in vectors:
        while r:
            p = r.bit_length() - 1
            if p in pivots:
                r ^= pivots[p]
            else:
                pivots[p] = r
                break
    def reduce(v):
        while v:
            p = v.bit_length() - 1
            if p not in pivots:
                return v
            v ^= pivots[p]
        return 0
    return pivots, reduce

_, reduce_trivial = make_reducer(trivial)

nondeg_bit = 1 << idx[(0, frozenset(), frozenset())]

# find a null basis vector outside the twist-family span
extra = None
for v in null_basis:
    if reduce_trivial(v):
        extra = v
        break
assert extra is not None
s = extra if (extra & nondeg_bit) else (extra ^ can_vec)
assert s & nondeg_bit, "should be non-degenerate now"
assert reduce_trivial(s), "should be outside twist family"

# unpack s into a construction dict
D = {}
for k, (n, V, W) in enumerate(unknowns):
    if s >> k & 1:
        i = n - len(V) - len(W)
        D.setdefault((i, n), set()).add((V, W))

def Ds(i, n):
    return D.get((i, n), set())

print(f"non-degenerate irreducible solution outside twist family (level <= {M}):")
print(f"  deviation from nearest twist? components differing from can and Tcan:")
for (i, n) in sorted(D.keys(), key=lambda t: (t[1], t[0])):
    ch = Ds(i, n)
    c = canonical(i, n)
    tags = []
    if ch != c and ch != T(c):
        tags.append("DIFFERS from can & Tcan")
    if not is_free(ch, i, n):
        tags.append("NOT FREE")
    assert is_irreducible(ch)
    if tags:
        print(f"  (i={i}, n={n}): {', '.join(tags)}; {len(ch)} terms")

# show the smallest deviating component explicitly
for (i, n) in sorted(D.keys(), key=lambda t: (t[1], t[0])):
    ch = Ds(i, n)
    c = canonical(i, n)
    if ch != c and ch != T(c):
        print(f"\n  explicit Delta_{i}[{n}] =")
        for V, W in sorted(ch, key=lambda p: (sorted(p[0]), sorted(p[1]))):
            print(f"    {tuple(sorted(V))} ot {tuple(sorted(W))}")
        print(f"  canonical_{i}[{n}] =")
        for V, W in sorted(c, key=lambda p: (sorted(p[0]), sorted(p[1]))):
            print(f"    {tuple(sorted(V))} ot {tuple(sorted(W))}")
        diff = ch ^ c
        print(f"  difference (Delta + canonical) =")
        for V, W in sorted(diff, key=lambda p: (sorted(p[0]), sorted(p[1]))):
            print(f"    {tuple(sorted(V))} ot {tuple(sorted(W))}")
        break
