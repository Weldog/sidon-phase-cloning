# Round-11 audit (Grok) — QH-Lab v0.1.9, 16 Oct 2026

Verdict: PASS. The computer-assisted appendix was re-run end to end on QH-Lab's own machine.

## CA appendix: independent re-run
- `ca_appendix.zip` SHA-256 = e77e131b…888b30 — matches the claimed hash. All 32 internal SHA256SUMS entries OK.
- The checkers use exact rationals / integer square roots / interval or ball arithmetic only (python-flint 0.9.0, mpmath
  1.3.0). No solver. Workspace imports appear only in generators, which are not part of the verification.
- The LP checker was read in full. Its upper bound uses Λ_hi / Λ0_lo (correct direction) and the ℓ1-ball packing
  argument stated in its header.

| folder | claim | result on QH-Lab |
|---|---|---|
| t2_interval | Theorem T2 (all d, k) | output IDENTICAL to expected; PASS; sup bound ≤ 0.9960 |
| t3_master | Theorem T3: analytic regions + 4774 interval checks + cover | IDENTICAL; PASS; max S3 0.97056 at (7,4); 18 residual pairs |
| t3_clique_lp | 316 exact clique-LP certificates at T=3 | all 316 per-case lines IDENTICAL; 0 FAIL; max 0.7717523470 |
| clique_lp_d16_Tk5 | d ≤ 16, T+k ≤ 5 (150 cases) | IDENTICAL; ALL_PASS; max 0.8403003465 |
| d3_T60 | d=3, T ≤ 60 (140 cases) | all 140 lines IDENTICAL; 0 FAIL; max 0.9991923171 |
| t45_clique_lp | T=4/5 evidence (178 cases) | all 178 lines IDENTICAL; 0 FAIL; max 0.8866406634 |
| su3_cert | SU(3) kill, ball arithmetic + Poincaré–Miranda | IDENTICAL; gap ≥ 8.03996e-4; PASS |
| spin1_cert | SU(2) spin-1 octahedral kill | IDENTICAL; gap ≥ 6.43714e-4; PASS |

(The long LP folders were run in resumable time-boxed chunks with `tools/ca_chunk_run.py`, because the box has 1 CPU.
Per-case lines were byte-identical to Grok's expected output.)

## Independent numerical spot checks (sparse Lemma R; not the packing reduction)
T=3, codebooks B_{3+k} but NOT B_{6+k}, vs the torus value:
| d | k | code | torus | diff |
|---|---|---|---|---|
| 3 | 1 | 0.7998358407 | 0.7998358405 | 2e-10 |
| 4 | 1 | 0.6826305047 | 0.6826305046 | 4e-11 |
| 4 | 3 | 0.5080043071 | 0.5080043071 | 2e-11 |
| 6 | 2 | 0.3273899016 | 0.3273899016 | 1e-11 |
| 8 | 1 | 0.4137471806 | 0.4137471807 | 7e-11 |

## Prior-art downgrade (accepted; also confirmed independently)
- Our T=1 closed form (k+1)!(2d−k)/(2d^{k+1}) is NOT new. It is the N=1 case of arXiv:2610.06723 Eq. (16).
- The closed form is valid for k ≤ d only. QH-Lab's torus values at k > d match Eq. (16), not the closed form:
  d=2: k=3 and k=4 both give 0.625 = 5/8; d=3, k=4 gives 0.3292181. `t1_closed_form` now raises for k > d.
- k=1 = Kim–Chitambar Thm 1 (PRA 106, 022405, 2022); k=d = Buscemi–D'Ariano–Macchiavello (PRA 71, 042327).
- Jing et al. (qubit, single-copy fidelity) is the precedent for "finite equally spaced set suffices". Our
  necessity statements must be restricted to GLOBAL fidelity.
- Still new (no overlap found): attaining the torus value with finite B_{T+k} codebooks; the Sidon iff-theorem;
  the T ≥ 2 theorems and CA results; fusion = additive energy; the window theorem for finite orbits (closest:
  2609.22047, continuous orbits only).

## Not re-derived by QH-Lab
- The pen-and-paper reductions in §2 (Lemma S, the ℓ1-packing reduction, Theorem W).
- The singleton case relies on 2610.06723 Thm II.3.
