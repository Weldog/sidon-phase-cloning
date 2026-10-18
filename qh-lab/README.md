# QH-Lab — Quantum Hash research harness (B1 core)

Reproducible numerics for QDH-1 (phase hash) and QDH-2 (design codebooks). Security model: Tier B
(unbounded adversary, standard QM), G4 copy limits. No hardware calls in B1.

## Run (Windows)
    cd qh-lab
    python -m venv .venv
    .venv\Scripts\activate
    pip install -r requirements.txt
    python -m pytest -q tests        # 49 acceptance tests, ~1-2 min
    python run_b1.py                 # writes results.json, ~30 s
    python r4_verify.py              # Round-4 audit, writes r4_verify.json, ~1 min
    python r5_verify.py              # Round-5 audit, writes r5_verify.json, ~1 min
    python bqsm_mask_attack.py       # QDH-BQSM no-go: zero-memory attack on shift-masked broadcast
    python r6_verify.py              # Round-6 audit: Theorem F, Conjecture A non-coherent tests (~5 min)
    python r7_verify.py              # Round-7 audit: spin-3/2 GHZ-like genuine test + controls
    python r8_verify.py              # Round-8 audit: SU(3)/Sym^2 Hessian counterexample (exact Schur-Weyl twirl)
    python r9_verify.py              # Round-9 audit: spin-1 octahedral kill + torus minimal-q spot checks
    python t2_spot.py                # Round-10 audit: T=2 Torus Theorem spot checks up to d=20 (sparse Lemma R)
    python r10_spin52.py             # Round-10 audit: spin-5/2 binary-icosahedral window prediction (fail)
    python tools/ca_chunk_run.py <ca_appendix/folder> 250   # resumable runner for Grok's exact LP checkers (1-CPU boxes)
    python tools/neg_table_check.py  # reproduces all 75 values of the paper's negative-results table
    python tools/coverage_check.py <path/to/ca_appendix>   # region-coverage audit (SHA-256 c42d11f6...)

## Modules
| file | contents |
|---|---|
| `qhlab/phase_hash.py` | QDH-1 codebook, FFT certificate of delta (O(q log q)), brute-force check, Welch/Levenstein floors |
| `qhlab/designs.py` | D_T, Werner cloning value W, Haar moment operator, design defect eta_{T,k}, frame potential |
| `qhlab/attacks.py` | optimal preimage SDP, PGM via Gram matrix, optimal counterfeiting SDP (Choi), binding bound |
| `qhlab/stabilizer.py` | n-qubit stabilizer states (exact 3-design) by Clifford BFS |
| `qhlab/clifford.py` | Clifford group n=1,2 (24, 11520), orbit codebooks, fiducial search on collision bias |
| `qhlab/theorems.py` | Round-4: Theorem A exact formula, Theorem B certificate C_{T,k} (replaces eta), Witting-40, Sidon test |
| `qhlab/lemma_r.py` | Round-5: independent Lemma R reduced dual, torus K, T=1 closed form, Theorem-1 Kraus channel, Bose–Chowla sets |
| `qhlab/fib.py` | Fibonacci braid generators — benchmark only (hash family KILLED, Round 3) |
| `qhlab/delaunay.py` | d=2 closest-pair-is-Delaunay-edge unit test only (certifier KILLED as an advantage) |

## Acceptance tests (tests/test_b1.py)
Fibonacci unitarity / Yang–Baxter / far commutation / sigma^10=1 (d=2,5,13,34); FFT delta == brute force;
stabilizer counts 6/60 and exact 3-design, not 4; SIC 2-design not 3; Haar moment trace; Round-3 numbers
reproduced; stabilizer counterfeiting == Werner; Delaunay closest-pair lemma (100 trials).
