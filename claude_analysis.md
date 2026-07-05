Let me read the key sections to understand what changed mathematically.

Now let me check the examples and legacy sections for additional context.

Here is a clear account of the mathematical evolution since the PEMS submission.

---

## Mathematical changes since `c6c0501`

### 1. The axiom "non-degenerate" was replaced by "non-zero"

In the submitted version the main theorem used three axioms: **non-degenerate** (forces $x \cup_0 x \neq 0$ for 0-simplices), **irreducible**, **free**. After the submission you replaced "non-degenerate" with **non-zero** (meaning "not identically zero"), which is strictly weaker. The theorem statement was also extended from *simplicial* to **semi-simplicial** cup-$i$ constructions — a natural generalization since semi-simplicial sets lack degeneracy maps.

### 2. A counter-example was found (branch `better`, commit `4d25956`)

You found a cup-$i$ construction that is **irreducible but not free**. This was constructed explicitly in `sec/examples.tex` by what you now call the "adding partition boundaries" move: starting from the canonical construction $\canonical$, you perturb $\canonical_0[2]$ by adding $\bd\zeta_{\{1\}}$. The result is irreducible but violates freeness. This confirmed that **freeness is genuinely independent** of the other axioms, and motivated the later attempt (branch `no_freeness`) to understand whether freeness could be removed from the *theorem statement* at all.

### 3. The proof was substantially restructured around "partition chains"

The old proof worked by a direct induction, using the kernel-of-codegeneracy condition and the functions $\xi: \mathcal{P}^n_{n-i} \to \mathbb{F}_2$ to nail down $\triangle_i[n]$ to either $\Delta_i[n]$ or $T\Delta_i[n]$. Two key lemmas ("first nail", "second nail") drove that argument.

The new proof introduces a new piece of intermediate structure: **partition chains**. For $U \subseteq \{0,\dots,n\}$, the partition chain $\zeta_U$ is the sum of all basis elements $V \otimes W$ in $\mathcal{P}(\Delta^n)^{\otimes 2}$ with $V \cup W = U$ and $V \cap W = \emptyset$. The key new lemma (`l:partition_chains`) is:

> *An irreducible chain $\zeta$ satisfies $(\pi_{\mathrm{red}} \circ \partial)\zeta = 0$ if and only if $\zeta$ is a sum of partition chains.*

The proof then proceeds in cleaner case-by-case stages:
- **$i \geq n$ case**: concluded quickly from freeness and the structure of $\mathcal{P}(\Delta^n)^{\otimes 2}$.
- **$i = n-1$ case**: the new argument writes $\triangle_{n-1}[n]$ as $\canonical_{n-1}[n]$ plus a correction that must be a sum of partition chains of singletons, and uses a boundary argument to show that correction must be zero or all of $(1+T)\canonical_{n-1}[n]$.
- **Induction step** (new $k = n-i > 1$): uses the partition chain lemma to show that any deviation from $\canonical_i[n]$ would require a partition chain of a set with more than one element, but the number of such partitions exceeds 2, giving a contradiction.

### 4. The old proof material was preserved in `sec/legacy.tex`

`sec/legacy.tex` grew from 635 to 914 lines. The old "first nail / second nail" lemmas are still there, flagged with `\anibal{Below the old}`. They have not been deleted; they are archived in legacy.

### 5. There are unresolved issues flagged in the current HEAD

Two `\anibal{...}` comments mark open problems:

- In `sec/proof.tex`: `\anibal{This proof doesn't use freeness}` — inside the $i = n-1$ case, suggesting freeness may be redundant there (which would be important if you want to weaken the axioms).
- In `sec/examples.tex`: `\anibal{Missing the $i_0$, $n_0+1$ case, for which it fails it seems}` — the "adding partition boundaries" move is stated without a proof covering all cases.

---

**In short**: the post-submission work is a genuine mathematical development, not cosmetic revision. The counter-example is real and the proof strategy was changed to use partition chains. But the current HEAD is explicitly unfinished work-in-progress, with acknowledged gaps. If the PEMS submission is still under review and you want to revert to a clean state, `c6c0501` is stable; if you want to continue with this new approach, the two `\anibal{}` flags are the live issues.