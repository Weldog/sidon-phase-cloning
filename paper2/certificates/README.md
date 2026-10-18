# Exact certificates for "Universal Cloning Is Optimal Exactly for Designs"

`exact_cert.py` certifies, in exact Gaussian-rational arithmetic, the optimal average global fidelity P_E(1,k) of T=1
cloning for finite ensembles whose amplitudes lie in Z[i]/2^(r/2):
- DUAL: Y = lam*I; positive semidefiniteness of I⊗Y − X is checked by an exact LDL factorization (real embedding).
- PRIMAL: an exact projector onto the lam-eigenspace, rescaled to a Choi operator; Tr_out J = I is checked exactly.
Primal value = dual value = d*lam, so the value is certified.

`h2_stabilizer.py` generates the state sets (stabilizer states by Clifford-orbit search from |0…0⟩, with no external
data). It was written independently of `exact_cert.py`.

Run (Python ≥ 3.10, numpy):  python3 exact_cert.py
Expected (see small_certs.json, stab14_cert.json):
  octahedron (6)            1→2 = 2/3, 1→3 = 1/2, 1→4 = 5/12   (universal 2/3, 1/2, 2/5)
  d=4 complete MUBs (20)    1→2 = 2/5, 1→3 = 1/4               (universal 2/5, 1/5)
  2-qubit stabilizer (60)   1→4 = 1/8                          (universal 4/35)
Each line ends "CERTIFIED=True". The 2-qubit stabilizer 1→4 case takes a few minutes.

Provenance: exact_cert.py was written by Grok Bot (xAI) and re-run independently by Claude (Anthropic) with the
independently written h2_stabilizer.py; both runs give identical values.
