# Paper 3 — code and certificates ("Exact Optimal Cloning of Stabilizer States")

certificates_grok/  (written by Grok Bot; re-run by Claude and on the author's PC with identical output)
- lagr.py            enumeration of the stochastic Lagrangian relations Sigma(q_{a,b}) and their composition
- alg.py             orbit-quotient counts of the left-multiplication matrix M (one change by Claude: the number of
                     worker processes is read from the PROCS environment variable, default 2, so t = 7 fits in 3 GB)
- cert_alg.py        exact symbolic certificates, valid for every n (D = 2^n)
- alg_<a>_<b>.json   counts for the 15 signatures with 3 <= t <= 7, a >= b, other than (3,3)
- cert_small.log, cert_t6.log, cert_t7.log   certificate output (Claude's run)
- stream_check.py    exact streaming check of Lemmas 7 and 8 over all 9,845,550 relations at t = 8
- verify_fixes.py, qidentity.py, consistency.py, predict_r26.py, orbit_structure_import.py   supporting checks
  (some contain absolute paths from the author's workspace; adjust before running)

independent_claude/  (written by Claude, sharing no code with certificates_grok/)
- thmC_claude.py     own enumeration of Sigma from the definition, own composition and matrix M: checks every column
                     sum (Lemma 7), every count (Lemma 8) and rho(M) = F_{t-w} at n = 1, 2, 3 for all 8 signatures
                     with a >= b and t <= 5; checks Lemma 9 for t <= 30
- witt_claude.py     exhaustive check of Proposition 6 for all 54 signatures with t <= 9
- r1g3_audit.py      brute-force checks: the many-copy identity, the T = 0 formula, spectral equivalences
- h2_stabilizer.py   generation of all n-qubit stabilizer states (n <= 3), used for direct Gram-kernel diagonalization
