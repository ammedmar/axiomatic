# Branch archived (2026-07-05)

This branch attempted to remove the freeness axiom from the main theorem, i.e.,
to prove that every non-zero irreducible semi-simplicial cup-i construction is
isomorphic to the canonical one. The rewrite of the induction step in
`sec/proof.tex` was left unfinished.

**The goal is impossible.** The construction

    (1+T)canonical = { canonical_i[n] + T canonical_i[n] }

is a valid simplicial cup-i construction (the defining identities are linear
and closed under composition with T), and it is:

- non-zero: (1+T)canonical_{n-1}[n] = sum_u ({u} ot {} + {} ot {u}) != 0;
- irreducible: transposition preserves irreducibility of basis elements;
- NOT free: e.g. {1} ot {} and {} ot {1} are both summands at (i,n) = (0,1);
- not isomorphic to the canonical construction, since cup-i isomorphisms
  (per-i twists by T) preserve freeness.

Hence non-zero + irreducible do not imply uniqueness, and freeness is a genuine
axiom. This counterexample now lives in `sec/examples.tex` on `main`
(subsection "Removing free"), replacing the earlier candidate counterexample
`canonical_0[2] + bd zeta_{1}`, which fails the defining identity at
(i,n) = (0,3) — the case flagged as missing in the proof of the
"adding partition boundaries" move.

Note: the observation recorded in this branch that the i = n-1 case does not
use freeness internally is correct; it survives on `main` as a remark after
that proof. Freeness enters that argument only through the i >= n case, where
it is essential: (1+T)canonical_n[n] = 0.

Whether freeness becomes redundant when the non-zero axiom is strengthened to
non-degeneracy (x cup_0 x != 0 for 0-simplices, which excludes (1+T)canonical)
remains open; computations in `python/claude_checks/` on `main` found
non-degenerate irreducible non-free solutions of the defining identities
truncated to n <= 6, suggesting the answer is negative there too.
