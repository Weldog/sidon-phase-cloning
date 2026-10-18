"""Round-7 audit (independent QH-Lab code): spin-3/2 GHZ-like icosahedral orbit at (1,1), tuned so M_2 matches SU(2).
Haar over SO(3) by an exact Euler product rule (validated); orientation tuned by least squares on M_2(X)-M_2(H)."""
import json, numpy as np, cvxpy as cp
from scipy.optimize import least_squares
from scipy.linalg import expm
j = 1.5; d = 4; m = np.arange(j, -j - 1, -1)
Jz = np.diag(m).astype(complex)
Jp = np.zeros((d, d), complex)
for i in range(1, d): Jp[i - 1, i] = np.sqrt(j * (j + 1) - m[i] * (m[i] + 1))
Jx, Jy = (Jp + Jp.T.conj()) / 2, (Jp - Jp.T.conj()) / 2j
def D_axis(axis, ang):
    a = np.asarray(axis, float); a = a / np.linalg.norm(a); return expm(-1j * ang * (a[0] * Jx + a[1] * Jy + a[2] * Jz))
def D_euler(al, be, ga): return expm(-1j * al * Jz) @ expm(-1j * be * Jy) @ expm(-1j * ga * Jz)
# icosahedral rotations (SO(3)) -> axis/angle -> spin-3/2
phi = (1 + 5 ** 0.5) / 2
def rot(axis, ang):
    a = np.asarray(axis, float); a /= np.linalg.norm(a); K_ = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K_ + (1 - np.cos(ang)) * K_ @ K_
G = [np.eye(3)]; gens = [rot([0, 1, phi], 2 * np.pi / 5), rot([1, 1, 1], 2 * np.pi / 3)]; fr = [np.eye(3)]
while fr:
    nx = []
    for A in fr:
        for g in gens:
            B = g @ A
            if not any(np.allclose(B, C, atol=1e-9) for C in G): G.append(B); nx.append(B)
    fr = nx
assert len(G) == 60
def so3_to_axis_angle(Rm):
    ang = np.arccos(np.clip((np.trace(Rm) - 1) / 2, -1, 1))
    if ang < 1e-9: return np.array([0, 0, 1.0]), 0.0
    if abs(ang - np.pi) < 1e-7:
        w, V = np.linalg.eigh((Rm + np.eye(3)) / 2); return V[:, -1], np.pi
    return np.array([Rm[2, 1] - Rm[1, 2], Rm[0, 2] - Rm[2, 0], Rm[1, 0] - Rm[0, 1]]) / (2 * np.sin(ang)), ang
DG = [D_axis(*so3_to_axis_angle(Rm)) for Rm in G]
# exact SO(3) Euler product rule
na, nb = 24, 24
al = 2 * np.pi * np.arange(na) / na; zb, wb = np.polynomial.legendre.leggauss(nb); be = np.arccos(zb)
Uq = [D_euler(a, b, c) for a in al for b, wbb in zip(be, wb) for c in al]
Wq = np.array([wbb / 2 / na / na for a in al for b, wbb in zip(be, wb) for c in al])
# validate: Haar average of D^{3/2} (x) conj D^{3/2} = swap/d  -> twirl of |0><0| equals I/d
twirl = sum(w_ * U[:, 0:1] @ U[:, 0:1].conj().T for U, w_ in zip(Uq, Wq))
assert np.allclose(twirl, np.eye(d) / d, atol=1e-12)
def kp(v, p):
    o = np.array([1.0 + 0j])
    for _ in range(p): o = np.kron(o, v)
    return o
def M(states, w, mm):
    V = np.array([kp(s, mm) for s in states]); return (V * w[:, None]).T @ V.conj()
def Rmix(states, w, T, k):
    V = np.array([np.kron(kp(s.conj(), T), kp(s, T + k)) for s in states]); return (V * w[:, None]).T @ V.conj()
def P_opt(states, w, T, k):
    X = Rmix(states, w, T, k); din, dout = d ** T, d ** (T + k)
    Jv = cp.Variable((din * dout, din * dout), hermitian=True)
    p = cp.Problem(cp.Maximize(cp.real(cp.trace(Jv @ X))), [Jv >> 0, cp.partial_trace(Jv, [din, dout], axis=1) == np.eye(din)])
    p.solve(solver="SCS", eps=1e-10, max_iters=500000); return float(p.value)
def refinedB(states, w, T, k):
    V = np.array([kp(s, T) for s in states]); r = np.linalg.matrix_rank(V, tol=1e-9)
    Mm = M(states, w, T + k); return float(r * np.linalg.eigvalsh((Mm + Mm.conj().T) / 2)[-1])
def main():
    out = {}
    for name, fid in [("GHZ-like (|3/2>+|-3/2>)/rt2", np.array([1, 0, 0, 1], complex) / np.sqrt(2)),
                      ("|3/2,1/2>", np.array([0, 1, 0, 0], complex))]:
        H_states = np.array([U @ fid for U in Uq])
        MH = M(H_states, Wq, 2)
        def orbit(e3): return np.array([Dg @ D_euler(*e3) @ fid for Dg in DG])
        w60 = np.ones(60) / 60
        def resid(e3):
            Dm = M(orbit(e3), w60, 2) - MH; return np.concatenate([Dm.real.ravel(), Dm.imag.ravel()])
        best = None
        for start in [(0.1, 0.4, 0.2), (0.7, 1.1, 0.3), (1.3, 0.6, 2.0), (2.2, 1.9, 0.9)]:
            r = least_squares(resid, start, xtol=1e-15, ftol=1e-15, gtol=1e-15)
            if best is None or r.cost < best.cost: best = r
        tuned = orbit(best.x); untuned = orbit((0.3, 0.7, 1.1))
        PH = P_opt(H_states, Wq, 1, 1)
        row = dict(fiducial=name, P_H=PH, thmB_H=refinedB(H_states, Wq, 1, 1))
        for tag, X in [("tuned", tuned), ("untuned", untuned)]:
            row[tag] = dict(M2_diff=float(np.abs(M(X, w60, 2) - MH).max()),
                            mixed_diff=float(np.abs(Rmix(X, w60, 1, 1) - Rmix(H_states, Wq, 1, 1)).max()),
                            P_X=P_opt(X, w60, 1, 1), thmB_X=refinedB(X, w60, 1, 1))
            row[tag]["gap"] = row[tag]["P_X"] - PH
        out[name] = row; print(json.dumps(row, indent=1, default=float), flush=True)
    json.dump(out, open("r7_verify.json", "w"), indent=1, default=float)

if __name__ == "__main__":
    main()
