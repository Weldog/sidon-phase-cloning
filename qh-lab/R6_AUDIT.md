# Round-6 audit (Grok) — QH-Lab v0.1.4, 16 Oct 2026

## Verified
| claim | Grok | QH-Lab (independent) |
|---|---|---|
| Theorem F: optimal fusion = E_q(K)/d³ | 7 sets | SDP = formula: [0,1,3]/7 0.555556; [0,1,2]/7 0.703704; [0,1,5,17]/32 0.46875; [0,1,2,3]/32 0.6875 |
| Kill B (homomorphic unclonable hash) | killed | accepted; Theorem F is the keeper |
| Kill C (physical-law noisy storage) | killed in Tiers A/B | accepted (mass-density noise → spin/degenerate encodings evade; QEC handles the rest; noise parameters order-of-magnitude) |
| arXiv:2610.06723 Herasymenko–Lin–Verduyn Lunel–Grinko, "Approximate cloning of structured pure states" (5 Oct 2026) | prior art | exists |
| arXiv:1412.5135 Ziatdinov, "Quantum hashing. Group approach" | prior art | exists |
| arXiv:2609.22047 Ćwikliński–Studziński, cloning on orthogonal orbits | adjacent | exists — and it is exactly the real-state (O(d)-orbit) setting of our new tests |

## Correction: Grok's "genuine non-abelian" evidence was trivial
All of Grok's SU(2) passes use **coherent** fiducials. For those, a refined Theorem B
P_X ≤ dim(span{ψ_x^⊗T}) · λ_max(M_{T+k}(X))
already forces P_X = P_H: span = Sym^{2jT}(C²) and M_{T+k} = Π/(2j(T+k)+1) for a spherical 2j(T+k)-design,
so P_X ≤ (2jT+1)/(2j(T+k)+1) = P_H ≤ P_X. Checked: icosahedral 12-point (1,1) certificate = 0.6 = P_H.
Coherent-orbit Conjecture A is a corollary (Theorem B + generalised Werner), not evidence.

## First genuinely non-trivial non-abelian tests (QH-Lab)
Fiducial |1,0⟩ of spin-1 = real unit vectors in R³ under SO(3). Theorem B is loose here (shown), so equality is not forced.
P_H from exact product quadrature (degree ≥ 2(2T+k)+2).

| orbit (line-design strength) | (T,k) | (T+k)-moment match | mixed match | P_X | P_H | gap | Thm-B bound for H |
|---|---|---|---|---|---|---|---|
| icosa 6 axes (5) | (1,1) | yes | no | 0.600000 | 0.600000 | 2e-14 | 1.00 (loose) |
| icosa 30-line orbit at the l=6 zero (9) | (2,1) | yes | no | 0.766347 | 0.766347 | 8e-13 | 1.20 (loose) |
| icosa 30-line orbit (9) | (1,3) | yes | no | 0.428571 | 0.428571 | 1e-12 | 0.60 (loose) |
| icosa 6 axes (5) | (1,2) | **no** | no | 0.600000 | 0.600000 | −2e-13 | 0.60 (tight) |

- Conjecture A ("moment match ⇒ P_X = P_H"): **survives 3/3 non-trivial tests.**
- The "coincidence-order hierarchy is EXACTLY T+k ≤ ord(X)" claim: **FALSE.** Row 4 is not moment-matched yet
  equals P_H (here Theorem B is tight for H, so equality is forced anyway). Only the "if" direction survives.
- Evidence classes now: U(d) designs (proved), coherent orbits (proved via Thm B), torus B_{T+k} (numerical, R5),
  real/O(3) non-coherent orbits (numerical, this audit).

## Credit check
- ePrint 2026/1834 "Quantum Pessiland" lists Boyang Chen and Tomoyuki Morimae; Grok's added "Yamakawa" is not on
  that page → Unverified.
- Grok's Kazan correction (Quantum Reports 2025 = Ablayev, Khadiev, Vasiliev, Ziiatdinov) accepted; the brief had
  conflated it with Zinnatullin et al. 2025.
