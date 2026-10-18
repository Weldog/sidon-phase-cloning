# Sidon Sets, Additive Energy, and the Optimal Cloning of Finite Phase-State Codebooks

**Michael Weldon** · Independent Researcher · ORCID [0009-0001-6119-0560](https://orcid.org/0009-0001-6119-0560)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23428376.svg)](https://doi.org/10.5281/zenodo.23428376)

This repository contains the paper and everything needed to verify it. The archived, citable record is on Zenodo:
**doi:10.5281/zenodo.23428376** (all versions; v2 = 10.5281/zenodo.23444464). The Zenodo record is the version of record.

## What the paper shows
For quantum phase-state codebooks |ψ_x⟩ = d^{-1/2} Σ_j ω^{x k_j}|j⟩ (the states behind quantum fingerprinting and hashing),
how well an adversary can clone or fuse the states is governed by additive combinatorics of the key set K:
- 1→2 global cloning reaches the continuous optimum **iff K is a Sidon set**;
- B_{T+k} key sets attain the continuous (torus) optimum (proved for d = 2 and for T ≤ 3; conjectured in general);
- for T = 1, an exact characterisation for every k (B_{1+k} is sufficient but not necessary);
- optimal m-fold fusion fidelity equals the higher additive energy E_m(K)/d^{m+1};
- for finite orbits of non-abelian groups, a window theorem plus certified SU(2) and SU(3) counterexamples to moment matching.

## Contents
| Path | What it is |
|---|---|
| `paper/` | PDF, LaTeX source and bibliography of v2 |
| `ca_appendix/` | Computer-assisted certificates with standalone exact / interval / ball-arithmetic checkers (no optimisation solver). `SHA256SUMS` inside. Archive hash of the Zenodo copy `ca_appendix.zip`: e77e131b…888b30 |
| `coverage_check.py` | Region-coverage audit for the T = 2, 3 proofs. Usage: `python3 coverage_check.py ca_appendix` |
| `qh-lab/` | QH-Lab v1.1.0: independent numerical re-implementation, audit logs (R4–R11), 49 tests |
| `requirements.txt` | Exact library versions used for the published runs |

## Verify
```bash
pip install -r requirements.txt
cd ca_appendix && sha256sum -c SHA256SUMS && cat README.md   # then run each folder's checker
cd ../qh-lab && python -m pytest -q tests                      # expected: 49 passed
cd .. && python3 coverage_check.py ca_appendix                 # expected: "ALL_OK": true
```

## AI use
This work was carried out with substantial assistance from AI systems (Grok Bot, Claude, ChatGPT, Perplexity), as disclosed
in the paper's "Use of AI tools" section. Every computational claim is backed by an exact or interval certificate re-run
independently. The author takes full responsibility for the content.

## Licence
Code: MIT (see `LICENSE`). Paper: CC BY 4.0.
