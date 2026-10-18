"""Round-9 audit (independent QH-Lab code).
(1) SU(2)/SO(3) kill: spin 1 = vector rep on C^3 (Cartesian); G = 24 octahedral rotations (signed permutations, det +1);
    H = SO(3). Exact H-moments by an Euler product rule (own node counts, validated). P = 3*lmax(R) (Lemma 1) and a
    symmetry-free SDP cross-check. Grok's rational fiducial + own independently tuned fiducials.
(2) Torus spot checks with QH-Lab Lemma R: d=2 and d=3, B_{T+k} codebooks at minimal q (NOT B_{2T+k}) vs torus value."""
import itertools, json, numpy as np, cvxpy as cp
from fractions import Fraction as Fr
from scipy.optimize import least_squares
from qhlab import lemma_r as LR
# ---- octahedral rotations ----
O = []
for p in itertools.permutations(range(3)):
    for s in itertools.product((1, -1), repeat=3):
        Mx = np.zeros((3, 3))
        for i in range(3): Mx[i, p[i]] = s[i]
        if np.isclose(np.linalg.det(Mx), 1): O.append(Mx)
assert len(O) == 24
# ---- exact SO(3) product rule ----
def Rz(a): c, s = np.cos(a), np.sin(a); return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
def Ry(b): c, s = np.cos(b), np.sin(b); return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
na, nb = 16, 16
al = 2 * np.pi * np.arange(na) / na; zb, wb = np.polynomial.legendre.leggauss(nb)
Q = [Rz(a) @ Ry(np.arccos(z)) @ Rz(c) for a in al for z in zb for c in al]
Wq = np.array([w / 2 / na / na for a in al for w in wb for c in al])
# validation: average of R (x) R equals projector onto invariant vector |Phi><Phi| (1-dim), i.e. E[R_ij R_kl] = d_ik d_jl /3
E2 = sum(w * np.kron(R, R) for R, w in zip(Q, Wq))
assert np.allclose(E2, np.outer(np.eye(3).ravel(), np.eye(3).ravel()) / 3, atol=1e-13)
def kp(v, p):
    o = np.array([1.0 + 0j])
    for _ in range(p): o = np.kron(o, v)
    return o
def M2(states, w): V = np.array([kp(s, 2) for s in states]); return (V * w[:, None]).T @ V.conj()
def Rm(states, w): V = np.array([np.kron(s.conj(), kp(s, 2)) for s in states]); return (V * w[:, None]).T @ V.conj()
def lmax(A): return float(np.linalg.eigvalsh((A + A.conj().T) / 2)[-1])
def evaluate(f):
    X = np.array([g @ f for g in O]); H = np.array([R @ f for R in Q]); wX = np.ones(24) / 24
    RX, RH = Rm(X, wX), Rm(H, Wq)
    return dict(M2_diff=float(np.abs(M2(X, wX) - M2(H, Wq)).max()), mixed_diff=float(np.abs(RX - RH).max()),
                P_X=3 * lmax(RX), P_H=3 * lmax(RH), thmB=3 * lmax(M2(X, wX))), RX, RH
def sdp(R):
    J = cp.Variable((27, 27), hermitian=True)
    p = cp.Problem(cp.Maximize(cp.real(cp.trace(J @ R))), [J >> 0, cp.partial_trace(J, [3, 9], axis=1) == np.eye(3)])
    p.solve(solver="CLARABEL"); return float(p.value)
def main():
    out = {}
    fq = [(Fr(107192802089269, 200000000000000), Fr(-319249234016901, 500000000000000)),
          (Fr(402234308639641, 1000000000000000), Fr(-111070735693333, 500000000000000)),
          (Fr(-298671263020807, 1000000000000000), Fr(-34345633323871, 500000000000000))]
    f = np.array([float(a) + 1j * float(b) for a, b in fq]); f /= np.linalg.norm(f)
    r, RX, RH = evaluate(f); r["sdp_P_X"], r["sdp_P_H"] = sdp(RX), sdp(RH); r["gap"] = r["P_X"] - r["P_H"]
    out["grok_fiducial"] = r; print("Grok fiducial:", r, flush=True)
    own = []
    for seed in range(4):
        x0 = np.random.default_rng(500 + seed).normal(size=6)
        phi_of = lambda x: (lambda p: p / np.linalg.norm(p))(x[:3] + 1j * x[3:])
        def resid(x):
            X = np.array([g @ phi_of(x) for g in O]); H = np.array([R @ phi_of(x) for R in Q])
            Dm = M2(X, np.ones(24) / 24) - M2(H, Wq); return np.concatenate([Dm.real.ravel(), Dm.imag.ravel()])
        f2 = phi_of(least_squares(resid, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000).x)
        r2, _, _ = evaluate(f2); r2["gap"] = r2["P_X"] - r2["P_H"]; own.append(r2); print("own seed", seed, r2, flush=True)
    out["own_fiducials"] = own
    # ---- torus spot checks (Lemma R): minimal-q B_{T+k} codebooks vs torus value ----
    tor = []
    for d, T, k in [(2, 5, 1), (2, 8, 2), (2, 12, 1), (2, 15, 3), (3, 4, 1), (3, 6, 2), (3, 8, 1), (3, 3, 3)]:
        h = T + k
        K = [0, 1] if d == 2 else [0, 1, h + 1]
        q = h + 1 if d == 2 else h * (h + 1) + 1
        assert LR.is_Bh(K, q, h) and not LR.is_Bh(K, q, 2 * T + k)
        Kt, qt = LR.torus_K(d, T, k)
        tor.append(dict(d=d, T=T, k=k, K=K, q=q, P_code=LR.lemma_r(K, q, T, k), P_torus=LR.lemma_r(Kt, qt, T, k)))
        tor[-1]["diff"] = tor[-1]["P_code"] - tor[-1]["P_torus"]; print(tor[-1], flush=True)
    out["torus_spot"] = tor
    json.dump(out, open("r9_verify.json", "w"), indent=1, default=float)
if __name__ == "__main__":
    main()
