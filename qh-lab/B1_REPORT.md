# QH-Lab B1 report — 15 Oct 2026

**Status:** B1 core built and passing (16/16 tests). `run_b1.py` runtime 28 s. All numbers below are
from `results.json` (seed 20261015). Labels: Established / Derived (proof or numerics here) / Conjecture.

## A. Round-3 (Grok) reproduction — 7/7 match
| claim | Grok | QH-Lab |
|---|---|---|
| SIC d=2 preimage T=2 | 0.750 | 0.7500 |
| SIC d=2 preimage T=3 | 0.9714 | 0.9714 |
| Phase d=2 N=8 preimage T=2 | 0.3643 | 0.3643 |
| Phase d=2 N=8 counterfeit 1→2 | 0.750 | 0.7500 |
| Phase d=2 N=16 counterfeit 1→2 | 0.750 | 0.7500 |
| Phase d=2 N=16 counterfeit 1→3 | 0.750 | 0.7500 |
| Phase d=2 N=8 design defect η₁,₁ | 0.1405 | 0.1405 |

## B. QDH-2 stabilizer arm (new)
| n | d | N | L bits | δ² | 3-design | 4-design | preimage PGM = bound for | counterfeit vs Werner W |
|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 6 | 2.6 | 0.5 | yes | no | T ≤ 3 | 1→2, 1→3, 2→3 all = W |
| 2 | 4 | 60 | 5.9 | 0.5 | yes | no | T ≤ 3 | 1→2 = W (0.4); 1→3 = W (0.2) |
| 3 | 8 | 1080 | 10.1 | 0.5 | yes | no | T ≤ 3 | 1→2 = W (η=0 ⇒ theorem) |

- **Derived:** stabilizer codebooks remove QDH-1's counterfeiting excess at T=k=1 (exact, by η₁,₁=0).
- **Derived (numerics), Conjecture (general):** counterfeiting stays exactly at W even where η≠0
  (1→3 at n=1,2; 2→3 at n=1). The η upper bound is loose there. Conjecture S-W: for stabilizer
  codebooks, optimal T→T+k counterfeiting equals Werner's W for all n, T, k. Needs proof (Clifford
  commutant / Gross–Nezami–Walter structure).
- **Derived — the price:** collision bias is stuck at δ² = 1/2 for every n (distinct stabilizer states
  overlap at most 1/2). Binding per instance ≤ (1+2^{-k/2})/2 → needs many verification copies k.
  Index length is only L ≈ n²/2 bits, and preimage breaks at T* ≈ n/2 (consistent with stabilizer
  learnability from O(n) copies).

## C. Phase hash at d=4 (QDH-1 control)
| L | N | best δ² (200 K draws) | η₁,₁ | counterfeit 1→2 | W | Levenstein floor |
|---|---|---|---|---|---|---|
| 6 | 64 | 0.636 | 0.0345 | 0.4375 | 0.400 | 0.360 |
| 13 | 8192 | 0.983 | 0.0345 | 0.4375 | 0.400 | 0.400 |
- **Derived:** the flat-modulus counterfeiting excess (+0.0375 over W) is N-independent at d=4 too,
  extending Grok's d=2,3 finding.

## D. Clifford-orbit arm (new)
- **Derived:** any fiducial's orbit under the Clifford group is an exact 3-design (η₁,₁ = 0,
  counterfeit 1→2 = W at n=1).
- **Derived:** best fiducial found gives δ² = 0.862 (n=1, N=24) and 0.926 (n=2, N=11520), far above
  Levenstein (0.636, 0.400). Large orbits buy index length (L ≈ 13.5 bits at n=2) at the cost of
  collision resistance.

## Verdict
There is a three-way trade-off inside G4: **design order (counterfeiting) vs. collision bias δ vs.
index length L**. Phase codes win on δ and L but leak a constant counterfeit excess. Stabilizer codes
hit the cloning optimum but cap δ² at 1/2 and L at ~n²/2. Clifford orbits get L but lose δ. No arm
wins all three. No new resource — Proposition 1 stands.

## Open for Round 4
1. Prove or break Conjecture S-W.
2. Is there a codebook family with design order ≥ 3, δ² → 0, and L ≫ n² (e.g. orbits of a
   larger unitary 3-design, or approximate designs from random circuits)?
3. Does a product/concatenated codebook (stabilizer ⊗ phase, same index x) inherit W-level
   counterfeiting and phase-level δ?
