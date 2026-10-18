# Round-7 audit (Grok) — QH-Lab v0.1.5, 16 Oct 2026

Verdict: PASS. All numerics reproduced with independent QH-Lab code (exact SO(3) Euler product rule, validated;
orientation tuned by least squares on M_2(X) − M_2(H)).

| Grok row | case | Grok | QH-Lab |
|---|---|---|---|
| 7 | spin-3/2 GHZ-like, tuned, (1,1) | P_X = P_H = 0.442857143, Thm B 0.6, mixed diff 5.6e-3 | P_X = P_H = 0.442857143 (= 31/70), gap 1.6e-13; M_2 diff 2e-16; mixed diff 7.2e-3; Thm B 0.6 → **genuine pass** |
| 8 | spin-3/2 |3/2,1/2⟩, untuned control (0.3,0.7,1.1) | 0.419437755 vs 0.419047619, gap 3.9e-4 | 0.419437755 vs 0.419047619, gap 3.90e-4 (identical to 9 digits) |
| 6 | spin-3/2 |3/2,1/2⟩, tuned | mixed moment also matched → trivial | confirmed (mixed diff 6e-16) |
| extra | GHZ-like, untuned | — | gap 0.0202 → under-designed orbits fail strictly |

Prior art (IDs verified): Chiribella–Yang arXiv:1404.0990 "Optimal asymptotic cloning machines"; Kay arXiv:1407.4951
"Optimal Universal Quantum Cloning: Asymmetries and Fidelity Measures". Verdict PARTIALLY IMPLIED accepted:
coherent/highest-weight and U(d) cases are covered (refined Theorem B = Chiribella–Yang Eq. 12); non-coherent finite
orbits and finite torus subsets are not.

**Correction of our own audit:** ePrint 2026/1834 metadata lists Boyang Chen, Tomoyuki Morimae, Takashi Yamakawa.
Grok was right; our R6 credit note was wrong (based on a truncated search snippet).

Not independently re-checked: torus computer-assisted proof (36 sets; floating-point MILP), Hatanaka break
(arXiv:2409.19932), Eq. (16) cross-check against 2610.06723.

## Conjecture A — evidence ledger
| class | status |
|---|---|
| U(d) designs | proved (Werner / Theorem B) |
| coherent / highest-weight orbits | proved (refined Theorem B; Chiribella–Yang) — prior art |
| torus, B_{T+k} | computer-assisted for 36 (d,T,k), any q (float MILP); 660 random codebooks + R5 sweep numerical |
| real O(3) orbits, non-coherent |1,0⟩ | 3 genuine passes (R6 audit) |
| SU(2) spin-3/2 GHZ-like | 1 genuine pass (R7, reproduced) |
| SU(3) Sym²/adjoint | not run |
| counterexamples | none |
