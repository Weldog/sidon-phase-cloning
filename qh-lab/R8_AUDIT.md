# Round-8 audit (Grok) — QH-Lab v0.1.6, 16 Oct 2026

## Headline: Conjecture A is FALSE — confirmed independently
Independent construction: Hessian group = qutrit Clifford group mod phases (BFS from X, Z, F, diag(1,1,ω): order 216;
Σ|tr|²/|G| = 1.000 → irreducible on Sym²C³). Exact SU(3) moments by Schur–Weyl: twirl under U^{⊗n} = HS projection onto
span{P_σ : σ∈S_n} (pseudo-inverse Gram; n=4 for M_2, n=6 with partial transposes for R_H). Validation: invariance residual
2e-17 under random SU(3); tr R_H = 1; M_2(H) from the n=4 twirl equals the partial trace of the n=6 R_H to 6e-17.
Own fiducials (random starts, seeds 100–103), own least-squares tuning of M_2(X) = M_2(H).

| seed | M_2 diff | mixed diff | P_X | P_H | gap | refined Thm B |
|---|---|---|---|---|---|---|
| 0 | 3e-17 | 3.6e-3 | 0.3061037068 | 0.3060948064 | 8.9e-6 | 0.405 |
| 1 | 3e-17 | 3.2e-3 | 0.2894177299 | 0.2894068891 | 1.1e-5 | 0.336 |
| 2 | 3e-17 | 2.4e-3 | 0.2857753652 | 0.2857747813 | 5.8e-7 | 0.288 |
| 3 | 6e-17 | 4.0e-3 | 0.2857675396 | 0.2857673150 | 2.2e-7 | 0.292 |

- 4/4 moment-matched orbits fail strictly (gaps 2e-7 … 1e-5, i.e. 8–10 orders above numerical precision). Grok's
  headline fiducial gave 8.0e-4; its other two gave 9.7e-6 and 3.2e-5 — same phenomenon, fiducial-dependent size.
- Seed 0 cross-checked by a symmetry-free SDP (SCS): P_X = 0.306103706785, P_H = 0.306094806356 — equal to 6·λmax(R)
  to 2e-13, so Lemma 1 is not doing the work.
- O(4) spot check: 24-cell (1,1) = 0.4999999978, P_H via exact 600-cell 11-design rule = 0.4999999994 → pass confirmed.
- Lemma 1, Lemma 2 and the off-diagonal-coupling obstruction: proofs re-read, correct.
- Prior art arXiv:1108.4540 exists: Xiong et al., "General Quantum Key Distribution in Higher Dimension" (content claim not re-read). Torus exact-rational certificates (150 + 23 sets) and the T=1
  analytic proof: not re-run here.

## State of the claims
| claim | status |
|---|---|
| Conjecture A (general, non-abelian) | **FALSE** (SU(3)/Sym², Hessian orbits) |
| torus: B_{T+k} ⇒ P_X = P_torus | T=1 analytic, all k ≤ d; exact rational CA for 2≤d≤16, T+k≤5 (150 sets) + 23 larger-T sets |
| U(d) designs, coherent orbits | proved (prior art) |
| real O(3)/O(4), spin-3/2 GHZ, 24/600-cell | pass (numerical); mechanism (vanishing coupling) not explained |
| corrected non-abelian condition R_X Π = R_H Π | sufficient (T=1, G irreducible); necessity open |
