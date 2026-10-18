# Round-10 audit (Grok) — QH-Lab v0.1.8, 16 Oct 2026

Verdict: PASS on everything checked.

## Torus Theorem at T=2 (all d, k) — independent spot checks
QH-Lab Lemma R is now fully sparse (dict-convolved representation counts, sparse constraints; no length-q arrays),
so the torus value is computable for astronomically large q. Codebooks K = {(T+k+1)^j}, which are B_{2+k} but NOT
B_{4+k}, vs the torus value:

| d | k | code | torus | diff |
|---|---|---|---|---|
| 5 | 1 | 0.4878719491 | 0.4878719488 | 3e-10 |
| 5 | 3 | 0.2394399610 | 0.2394399610 | 3e-11 |
| 8 | 1 | 0.3300582507 | 0.3300582505 | 2e-10 |
| 8 | 2 | 0.1440072979 | 0.1440072979 | 7e-13 |
| 12 | 1 | 0.2297414123 | 0.2297414115 | 8e-10 |
| 12 | 2 | 0.0701146561 | 0.0701146560 | 1e-10 |
| 16 | 1 | 0.1760241330 | 0.1760241332 | 2e-10 |
| 20 | 1 | 0.1426244844 | 0.1426244844 | 1e-11 |

Torus values also match the R5 torus table (d=12 (2,2) 0.070115, d=16 (2,1) 0.176024, d=20 (2,1) 0.142624). The
interval-arithmetic proof (t2_proof.py, 1021 checks + 13 sharp cases) was not re-run here.

## Spin 5/2 under 2I — prediction replicated
Own implementation (SU(2)-lifted J-matrix rotations, exact 18³ Euler rule exact to L ≤ 17 ≥ 15 needed; 2I irreducible:
Σ|tr|²/60 = 1). 4/4 own tuned fiducials FAIL strictly, as the window rule predicts:

| seed | M_2 diff | mixed diff | P_X | P_H | gap | Thm B |
|---|---|---|---|---|---|---|
| 0 | 2e-16 | 2.4e-3 | 0.3553707680 | 0.3544893213 | 8.8e-4 | 0.402 |
| 1 | 1e-16 | 4.1e-3 | 0.3091298605 | 0.3089827116 | 1.5e-4 | 0.394 |
| 2 | 1e-15 | 3.6e-3 | 0.3114886669 | 0.3111828061 | 3.1e-4 | 0.324 |
| 3 | 2e-15 | 1.1e-3 | 0.5104318422 | 0.5103253227 | 1.1e-4 | 0.524 |

## Not re-derived here
Theorem W proof and the window formula ℓ ≤ 2min(a,b)+2T+k; O(4) window degrees; SU(3) ball-arithmetic re-certificate;
d=3 second implementation (T ≤ 60).

## Citations
arXiv:1111.4493 = Meagher–Purdy, "An Erdős–Ko–Rado theorem for multisets" — confirmed; cited for context only.
