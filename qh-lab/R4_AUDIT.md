# Round-4 audit (Grok) — QH-Lab v0.1.1, 16 Oct 2026

Verdict: PASS. Every checked number reproduced with QH-Lab's own code (not Grok's).

| claim | Grok | QH-Lab |
|---|---|---|
| S-W false: stab n=1 (1,3) | 0.416667 | SDP 0.416667, Thm A 0.416667 (W 0.4) |
| stab n=1 (2,2) / (3,1) | 0.625 / 0.833333 | SDP 0.625 / 0.833333 |
| stab n=2 (1,3),(2,2),(3,1) via Thm A | 0.125, 0.3125, 0.625 | identical |
| stab n=3 (1,3),(3,1) via Thm A | 0.027778, 0.416667 | identical |
| Thm B certificate = W iff T+k<=3 (n=2) | claimed | yes; (1,3) gives 0.2 > W 0.114 |
| Witting-40: N=40, delta^2=1/3 (Levenshtein-tight), 3- not 4-design, 1->2 = W | claimed | built independently (sign search, not Sp(4,3)); all confirmed |
| phase 1->2 = (2d-1)/d^2 "generic K" | d=4 | d=4,5,6 exact |
| additive coincidence K=[10,30,45,60] | 0.558 | 0.5578 |
| prior art arXiv:2609.20877 | Unverified | exists: "Stabilizer-Public-Key Authentication: Long Stabilizer Public Keys Resist Finite-Copy Forgery Attacks" |

Proof checks: Theorems A, B, C re-derived line by line — correct.

## New (QH-Lab, Derived)
"Generic K" = **Sidon K** (all pairwise sums k_i+k_j distinct mod 2^L).
- Sidon K: 1->2 counterfeiting = (2d-1)/d^2 exactly (d=4,5,6). Non-Sidon K leaks more (AP [0,1,2,3]: 0.745).
- P/W = 1 + (d-1)/(2d^2) -> 1. At d=2^n the excess is ~2^-(n+1) relative: Round-3's "phase family is weak"
  was a small-d artefact.
- Theorem C does not apply (phase codes are not Clifford-closed), and delta^2 -> 0 with eps-biased K.
- => Phase hash with K both Sidon and eps-biased is the only arm that is simultaneously near-W (1->2),
  delta^2 -> 0, and free in L. Only (T,k)=(1,1) is established; higher (T,k) is OPEN (B_h-set conjecture).
