# Computational companion

## `examples.ipynb`

Exploratory notebook based on the [`comch`](https://github.com/ammedmar/comch)
package (`pip install comch`). Contains the freeness/irreducibility checks for
the "transposing one element" example of the paper. Note: the "irreducible and
not free" cells document an earlier candidate counterexample that turned out to
violate the defining identity at (i, n) = (0, 3); the valid replacement is the
construction of `sec/examples.tex`, verified in `claude_checks/`.

## `claude_checks/`

Self-contained pure-Python scripts (no dependencies) verifying the paper's
results in the P-model of `sec/reformulation.tex`. Chains are sets of pairs of
frozensets (mod 2); the canonical construction passing all identities serves as
the sanity check for the model. Core library: `cupi.py`.

Main entry points, roughly in narrative order:

- `verify_axioms.py` — canonical construction sanity check; shows the two old
  candidate non-free examples fail at (0,3); validates the transposing example.
- `solution_space.py`, `extendable.py` — dimensions of the space of truncated
  irreducible constructions (2^{N+1} - N - 1) and surjectivity of restrictions.
- `simplicial_check.py` — imposing codegeneracy conditions changes no
  dimensions (numerical companion to "irreducible implies simplicial").
- `theorem_check.py` — exhaustive check that free + non-zero truncated
  solutions are exactly the twist family; verification of the partition chain
  lemma in all degrees for n <= 5.
- `kernel_check.py` — level-fresh deviations = augmented partition-chain
  cycles, dimension 2^n - 1.
- `single_column.py`, `stable_column.py`, `pair_extend.py`, `f_language.py`,
  `constants_space.py`, `c_functional.py` — the search narrowing down the
  single-column solutions to the constants triangle / Pascal transform
  description.
- `final_family.py` — verification of the closed-form family of
  `sec/examples.tex` (Theorem on the non-degenerate non-free construction)
  against the defining identities through n = 10.
