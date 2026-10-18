# Round-5 audit (Grok) — QH-Lab v0.1.2, 16 Oct 2026

Verdict: PASS. Lemma R was re-implemented from its statement alone (not Grok's core.py) and checked against
QH-Lab's unreduced Choi SDP.

| claim | Grok | QH-Lab (independent) |
|---|---|---|
| Lemma R = full SDP | 102 instances, ≤5.4e-8 | 8 cases across (1,1),(1,2),(2,1),(1,3),(3,1): ≤8.3e-8 |
| d=4 torus values, all six (T,k) | 0.4375 / 0.28125 / 0.577816 / 0.234375 / 0.428776 / 0.682631 | identical |
| T=1 closed form (k+1)!(2d−k)/(2d^{k+1}) | d≤20 | exact at d=4 (k=1,2,3), d=7 (k=2) |
| Theorem 1 primal channel gives (2d−1)/d² on ANY K | claimed | yes (Sidon, AP, d=5) |
| Sidon-only beats torus at T+k=3 | 25/25 cases | 3/3 random d=4 sets: (1,2) 0.292–0.326 vs 0.281; (2,1) 0.608–0.627 vs 0.578; (1,1) stays 0.4375; one confirmed by full SDP |
| Bose–Chowla p=7,h=3 | B3, δ²=0.525 | B3 yes, B4 no, δ²=0.525011; (1,2) = torus value; Theorem-2 certificate = formula |
| Bose–Chowla p=13,h=3 | δ²=0.296 | 0.296021 |
| Theorem 2 proof (B_m ⇒ E ψ^{⊗m} diagonal, λmax = m!/d^m) | — | re-derived, correct (needs d ≥ T+k) |

## Prior art Grok missed (verified)
- **Nagaj, Sattath, Brodutch, Unruh, "An adaptive attack on Wiesner's quantum money", QIC 16, 1048 (2016),
  arXiv:1404.1507.** Verification-oracle attacks (bomb-testing / protective measurement); fix = replace notes
  after validation. Round-5's oracle-leak result (§8.3) is this phenomenon for phase hashes → cite, do not
  claim as new.
- **Jing et al., "Minimal sets determining universal and phase-covariant quantum cloning", arXiv:1208.0746.**
  Qubit case: three equally spaced equatorial states already determine the phase-covariant cloner. Theorem 1 is
  the qudit generalisation with an exact additive-combinatorics (Sidon) characterisation → cite as the d=2
  precursor.

## Novelty after audit
Likely new (still Unverified until a fuller search): Theorem 1 (Sidon ⇔ phase-covariant optimum, qudits),
Lemma R, the active-sector criterion, Theorem 2 (short corollary), "ε-bias / Singer ETFs are worst for T+k≥3".
Known: states and ε-bias (Ablayev–Vasiliev), phase-covariant values for M=kd+N (Buscemi–D'Ariano–Macchiavello),
oracle leakage (Nagaj et al.), malleability (trivial).
