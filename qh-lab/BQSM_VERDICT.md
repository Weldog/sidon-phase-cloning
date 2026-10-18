# QDH-BQSM go/no-go — 16 Oct 2026

## 1. What Barhoush–Salvail actually use (read from source)
- "How to Sign Quantum Messages" (arXiv:2304.06325v5; Quantum 2026), Construction 4: the public key is a **program
  broadcast** of a random linear function P(x) = M·x (M a random m×n binary matrix). Signatures are two BQS
  one-time programs plus one-time quantum authenticated encryption of the message.
- Program broadcast and one-time programs come from "Powerful Primitives in the BQSM" (arXiv:2302.05724v6): built on
  **BQS oblivious transfer** (Damgård–Fehr–Salvail–Schaffner line), i.e. **BB84 conjugate-coding qubits**, entropic
  uncertainty relations, privacy amplification, and a new min-entropy splitting lemma. Honest one-time programs need
  only "prepare, send, and measure a single qubit in the BB84 basis".
- **No phase-hash / quantum-fingerprint states anywhere.** Novelty check: passed.

## 2. But QDH-BQSM fails on its own (Derived)
- Broadcasting many copies of a fixed hash state psi_x gives no bounded-storage protection. Measurement is not
  storage: an attacker with zero quantum memory measures each copy at once and keeps classical outcomes.
- Masking each copy with a fresh shift (psi_{x+r_i}, r_i revealed after t_end) makes every copy exactly I/d. It still
  fails. The shift-covariant POVM M_y = (d/q)|psi_y><psi_y| commutes with the mask, so the outcome y_i − r_i is a
  sample centred on x after the reveal.
- Simulation (`bqsm_mask_attack.py`): with **zero quantum memory**, x is recovered at 97.5% (d=8, L=10) and 100%
  (d=16, L=12) with T=64 copies.
- QDH's malleability, the same symmetry that makes it a valid hash family, is exactly what makes masking useless.
- The only repair is masking with **non-commuting (conjugate / Weyl–Heisenberg)** operators. Security then comes
  from the entropic uncertainty relation, which is the known high-dimensional BQSM, and the Sidon/B_h structure
  contributes nothing.

## Verdict
QDH-BQSM: **NO-GO.** The paper is the Sidon / B_h note alone.
