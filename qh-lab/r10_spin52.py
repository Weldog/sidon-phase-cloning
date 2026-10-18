"""Round-10 audit: spin 5/2 (d=6) under the binary icosahedral group, (T,k)=(1,1). Window rule predicts FAIL.
Own implementation: J-matrix rotations (SU(2) lift; sign irrelevant for states), exact SO(3) Euler rule (18^3, exact to
L<=17 >= 15 needed), own random fiducials tuned so M_2(X)=M_2(H). P = 6*lmax(R) (2I irreducible on V_5/2, checked)."""
import json, numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares
j = 2.5; d = 6; m = np.arange(j, -j - 1, -1)
Jz = np.diag(m).astype(complex); Jp = np.zeros((d, d), complex)
for i in range(1, d): Jp[i - 1, i] = np.sqrt(j * (j + 1) - m[i] * (m[i] + 1))
Jx, Jy = (Jp + Jp.conj().T) / 2, (Jp - Jp.conj().T) / 2j
def D_axis(a, ang): a = np.asarray(a, float) / np.linalg.norm(a); return expm(-1j * ang * (a[0] * Jx + a[1] * Jy + a[2] * Jz))
def D_euler(al, be, ga): return expm(-1j * al * Jz) @ expm(-1j * be * Jy) @ expm(-1j * ga * Jz)
phi = (1 + 5 ** 0.5) / 2
def rot(axis, ang):
    a = np.asarray(axis, float) / np.linalg.norm(axis); K_ = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K_ + (1 - np.cos(ang)) * K_ @ K_
G = [np.eye(3)]; fr = [np.eye(3)]; gens = [rot([0, 1, phi], 2 * np.pi / 5), rot([1, 1, 1], 2 * np.pi / 3)]
while fr:
    nx = []
    for A in fr:
        for g in gens:
            B = g @ A
            if not any(np.allclose(B, C, atol=1e-9) for C in G): G.append(B); nx.append(B)
    fr = nx
def axang(Rm):
    ang = np.arccos(np.clip((np.trace(Rm) - 1) / 2, -1, 1))
    if ang < 1e-9: return np.array([0, 0, 1.0]), 0.0
    if abs(ang - np.pi) < 1e-7: w, V = np.linalg.eigh((Rm + np.eye(3)) / 2); return V[:, -1], np.pi
    return np.array([Rm[2, 1] - Rm[1, 2], Rm[0, 2] - Rm[2, 0], Rm[1, 0] - Rm[0, 1]]) / (2 * np.sin(ang)), ang
DG = [D_axis(*axang(Rm)) for Rm in G]
irr = np.mean([abs(np.trace(U)) ** 2 for U in DG])
n = 18; al = 2 * np.pi * np.arange(n) / n; zb, wb = np.polynomial.legendre.leggauss(n); be = np.arccos(zb)
Ea, Eb = [expm(-1j * a * Jz) for a in al], [expm(-1j * b * Jy) for b in be]
Uq = np.array([A @ B @ C for A in Ea for B in Eb for C in Ea]); Wq = np.array([w / 2 / n / n for a in al for w in wb for c in al])
assert np.allclose(np.einsum('q,qi,qj->ij', Wq, Uq[:, :, 0], Uq[:, :, 0].conj()), np.eye(d) / d, atol=1e-13)
def kp(v, p):
    o = np.ones(1, complex)
    for _ in range(p): o = np.kron(o, v)
    return o
def M2(S, w): V = np.einsum('ni,nj->nij', S, S).reshape(len(S), -1); return (V * w[:, None]).T @ V.conj()
def Rm_(S, w): V = np.einsum('ni,nj,nk->nijk', S.conj(), S, S).reshape(len(S), -1); return (V * w[:, None]).T @ V.conj()
def lmax(A): return float(np.linalg.eigvalsh((A + A.conj().T) / 2)[-1])
w60 = np.ones(60) / 60
def main():
    out = {"irreducibility": float(irr), "rows": []}
    for seed in range(4):
        x0 = np.random.default_rng(900 + seed).normal(size=12)
        f_of = lambda x: (lambda p: p / np.linalg.norm(p))(x[:6] + 1j * x[6:])
        def resid(x):
            f = f_of(x); Dm = M2(np.array([U @ f for U in DG]), w60) - M2(Uq @ f, Wq)
            return np.concatenate([Dm.real.ravel(), Dm.imag.ravel()])
        f = f_of(least_squares(resid, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000).x)
        X, H = np.array([U @ f for U in DG]), Uq @ f
        RX, RH = Rm_(X, w60), Rm_(H, Wq)
        row = dict(seed=seed, M2_diff=float(np.abs(M2(X, w60) - M2(H, Wq)).max()), mixed_diff=float(np.abs(RX - RH).max()),
                   P_X=6 * lmax(RX), P_H=6 * lmax(RH), thmB=6 * lmax(M2(X, w60)))
        row["gap"] = row["P_X"] - row["P_H"]; out["rows"].append(row); print(row, flush=True)
    json.dump(out, open("r10_spin52.json", "w"), indent=1, default=float); print("irreducibility", irr)
if __name__ == "__main__":
    main()
