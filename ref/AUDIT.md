# Audit of the PEMS referee reply against the current manuscript

Date: 2026-07-06. The reply (`ref/referee.tex`, committed Jan 2024 in `c6c0501`)
answered a 14-item report on the PEMS submission. The manuscript has since been
substantially rewritten (semi-simplicial generalization, non-zero axiom,
partition-chain proof, independence section, simpliciality theorem), so this
audit checks each item against the CURRENT text, not the 2024 revision.

Legend: OK = fix present in current text; SUPERSEDED = the passage the referee
commented on no longer exists, and its replacement does not have the issue;
RESOLVED+ = went beyond what the referee asked.

| # | Referee point | Status | Where in current text |
|---|---|---|---|
| 1 | cup_i defined on cochains, applied to chains | OK | `statement.tex`, Axioms subsection: "Below we identify simplices with their associated chain and cochain basis elements", preceding the canonical definition |
| 2 | "that is natural" grammar | SUPERSEDED | sentence rewritten; current phrasing in `statement.tex` is grammatical |
| 3 | Lemma 14 "missing from V or W" ambiguity | SUPERSEDED | codegeneracy discussion rewritten; the sharp version is now in the proof of `t:irreducible_implies_simplicial` ("one of j, j+1 belongs to V and the other to W") and `eq:codegeneracies` was corrected to include the index shift |
| 4 | Lemma 17 "an element" of which set | SUPERSEDED | `ss:freeness_revisited` restated with explicit types (two-part statement) |
| 5 | Lambda(i,n) appears undefined | OK | `ss:axioms_revisited` introduces Lambda(i,n) explicitly |
| 6 | Lemma 17(1): how is naturality used | SUPERSEDED | the old lemma was replaced; the corresponding argument is now the i >= n case (`ss:i_geq_n`), rewritten with both transition cases |
| 7 | Lemma 19: is the partition unique | SUPERSEDED | the Lambda_1/Lambda_2 formulation was replaced by the xi-function formulation, where the choice is explicit |
| 8 | Reminder of canonical-vs-inspected notation | OK | statements in `proof.tex` carry `T^eps canonical` explicitly; `ss:proof` re-introduces the notation |
| 9 | Undeclared r in Section 4.7 | OK | no stray symbols found; `r` in the canonical definition is declared by the notation `U = {u_1 < ... < u_r}` |
| 10 | Lemma 8 or 9 citation | SUPERSEDED | old numbering gone after rewrite |
| 11 | Theorem 28 (Steenrod agreement) needs a proper proof | OK | `t:steenrod cup-i` in `others.tex` has a complete proof (non-zero, irreducible, free checks + the tie-break at [i+1]) |
| 12 | E_infty-operad definition should require freeness/projectivity | OK (no change needed, as argued in the reply) | `operads.tex`: definition kept as "same homology as Z, symmetric action free in each arity" |
| 13 | Cite Appendix 1 of the prop paper for the surjection operad | OK | `operads.tex`: "\cite[Appendix 1]{medina2020prop1}" |
| 14 | Suggestion: give examples showing the axioms (free, non-degenerate, irreducible) are independent | RESOLVED+ | previously marked NOT DONE. Now the entire Section `s:examples`: zero construction (non-zero), transposing example (irreducible), (1+T)canonical (free), plus `t:non_deg_non_free` establishing independence of freeness under the referee's own stronger non-degeneracy axiom, and `t:irreducible_implies_simplicial` as a by-product |

## Notes for the new cover letter / reply

- Item 14 is the headline: the referee's suggestion is now a full section with
  two new theorems, one answering the question for the stronger axiom system
  the referee's phrasing implicitly used (non-degenerate).
- Since the 2024 revision, the paper changed in ways the referee has not seen:
  (a) the theorem is now stated for semi-simplicial constructions;
  (b) the axiom "non-degenerate" was weakened to "non-zero" (the non-degenerate
  version is recovered as a corollary in `r:non_degenerate`);
  (c) the proof was restructured around partition chains;
  (d) new theorem: irreducible implies simplicial;
  (e) new section: independence of the axioms.
  A resubmission letter should list these as "changes beyond the report".
- "Other changes" from the 2024 reply (affiliation, referee acknowledgment):
  affiliation is current (Western University); `acknowledgment.tex` should be
  checked for the referee thanks before submitting.
