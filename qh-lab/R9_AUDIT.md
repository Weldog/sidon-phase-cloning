# Round-9 audit (Grok) — QH-Lab v0.1.7, 16 Oct 2026

## Headline: SU(2)/SO(3) kill — confirmed independently
Spin 1 as the vector rep on C³; G = 24 octahedral rotations (signed permutations, det +1); H = SO(3) via an own Euler
product rule (16×16×16, validated: E[R⊗R] = |Φ⟩⟨Φ|/3 to 1e-13). P = 3·λmax(R) (Lemma 1) plus a symmetry-free SDP.

| fiducial | M_2 diff | mixed diff | P_X | P_H | gap | Thm B |
|---|---|---|---|---|---|---|
| Grok's rational fiducial | 1.7e-16 | 1.4e-2 | 0.5165207341 (SDP 0.5165207339) | 0.5158770197 (SDP 0.5158770198) | **6.437e-4** | 0.682 |
| own seed 0 | 7e-17 | 1.5e-2 | 0.5002861079 | 0.5002715739 | 1.5e-5 | 0.505 |
| own seed 1 | 1e-16 | 4.6e-3 | 0.5700557193 | 0.5697816339 | 2.7e-4 | 0.581 |
| own seed 2 | 1e-16 | 9.5e-3 | 0.5796874863 | 0.5785912252 | 1.1e-3 | 0.934 |
| own seed 3 | 5e-16 | 7.1e-3 | 0.5858936294 | 0.5858011650 | 9.2e-5 | 0.957 |

5/5 moment-matched octahedral orbits fail strictly. Moment matching is NOT sufficient for SU(2) either. The
icosahedral passes (R6/R7) are consistent with Grok's window rule (first icosahedral invariant L=6 > window 4j=4).

## Torus: independent minimal-q spot checks (strongest regime: B_{T+k} but NOT B_{2T+k})
| d | T | k | K | q | code | torus | diff |
|---|---|---|---|---|---|---|---|
| 2 | 5 | 1 | {0,1} | 7 | 0.975217292 | 0.975217292 | 5e-10 |
| 2 | 8 | 2 | {0,1} | 11 | 0.991835602 | 0.991835601 | 6e-10 |
| 2 | 12 | 1 | {0,1} | 14 | 0.996649204 | 0.996649203 | 3e-10 |
| 2 | 15 | 3 | {0,1} | 19 | 0.992489903 | 0.992489898 | 5e-9 |
| 3 | 4 | 1 | {0,1,6} | 31 | 0.863874291 | 0.863874291 | 1e-11 |
| 3 | 6 | 2 | {0,1,9} | 73 | 0.907006826 | 0.907006824 | 1e-9 |
| 3 | 8 | 1 | {0,1,10} | 91 | 0.966372133 | 0.966372130 | 3e-9 |
| 3 | 3 | 3 | {0,1,7} | 43 | 0.733645989 | 0.733645988 | 5e-10 |
All agree to solver tolerance. Exact-rational certificates and the d=2 analytic proof were not re-run here.

## QH-Lab fix (own error)
`qhlab/lemma_r.py` was ill-conditioned for large T (unscaled counts up to d^(2T+k)); it returned garbage (>1) at
T ≥ 8. Fixed with column scaling + log-domain normalisation. All earlier QH-Lab results were at small T and are
unaffected (all R5 tests still pass).

## Citations (verified)
- arXiv:1404.1507 first author is Brodutch (Brodutch, Nagaj, Sattath, Unruh) — Grok's correction accepted.
- Gavinsky–Ito "Quantum Fingerprints that Keep Secrets" = arXiv:1010.5342 — confirmed.
