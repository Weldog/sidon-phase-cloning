# QHash R11 computer-assisted (CA) appendix

Every CA claim in qhash_round11.md has its own folder here. A folder holds:
- certificate data (JSON);
- a minimal standalone checker that uses exact rationals (Python `fractions`), integer square roots, or interval/ball arithmetic (mpmath `iv`, python-flint `arb`). No LP or SDP solver is used anywhere;
- `expected_output.txt`.

`SHA256SUMS` covers every file. Verify it with `sha256sum -c SHA256SUMS`.

| folder | claim | run (inside folder) | runtime |
|---|---|---|---|
| t2_interval | Theorem T2 (T=2, all d,k): constants C1-C3, 1021 finite interval checks C4, 13 Lagrangian residuals C5 | `python3 check_t2.py` | ~5 s |
| t3_master | Theorem T3 (T=3, all d>=3,k>=1): analytic-region constants, m0 table, 4774 interval checks of S3, and the cover check that every residual pair is LP-certified | `python3 t3_bounds.py 100 12 && python3 check_t3_cover.py` | ~3 s |
| t3_clique_lp | exact typed clique-LP certificates at T=3: d=2..12 with k<=20, and d=13..20 with k<=12 (316 cases, contains all 18 T3 residual pairs) | `python3 check_torus_lp.py cert.json` | ~3 min |
| clique_lp_d16_Tk5 | torus packing inequality for d<=16, T+k<=5 (150 cases) | `python3 check_torus_lp.py cert.json` | ~3 s |
| d3_T60 | d=3: T<=40 with k<=3, plus T=41..60 with k=1 (140 cases) | `python3 check_torus_lp.py cert.json` | ~1 min |
| t45_clique_lp | evidence for the T>=4 conjecture: T=4 (d=3..10 with k<=12; d=11..14 with k<=10) and T=5 (d=3..9 with k<=6), 178 cases | `python3 check_torus_lp.py cert.json` | ~2 min |
| su3_cert | SU(3)/Sym^2 kill: gap P_X-P_H >= 8.0399e-4 (ball arithmetic, Poincare-Miranda) | `python3 check_su3.py` | ~4 s |
| spin1_cert | SU(2) spin-1 / octahedral kill: gap >= 6.437e-4 (ball arithmetic, IVT) | `python3 check_spin1.py` | ~1 s |
| common | the shared clique-LP checker (identical copies sit in each LP folder) | - | - |
| generators | scripts that PRODUCED the certificates (HiGHS LP, float solvers). They are not part of the verification and are included for reproducibility only | - | - |

Requirements: Python >= 3.10, mpmath (t2_interval, t3_master), python-flint 0.9.0 (su3_cert, spin1_cert). The LP checkers need only the standard library.
