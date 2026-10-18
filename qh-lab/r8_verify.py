"""Round-8 audit: independent check of the SU(3)/Sym^2 Hessian counterexample to Conjecture A.
Exact SU(3) Haar moments via Schur-Weyl: twirl under U^{(x)n} = HS projection onto span{P_sigma, sigma in S_n}
(pseudo-inverse Gram, valid for n > d). R_H uses partial transposes for the conjugate factors. Own tuning, own seeds."""
import itertools, json, numpy as np
from scipy.optimize import least_squares
w3 = np.exp(2j * np.pi / 3)
# ---- Sym^2(C^3) isometry V: C^6 -> C^9 (real) ----
pairs = [(i, j) for i in range(3) for j in range(i, 3)]
V = np.zeros((9, 6))
for c, (i, j) in enumerate(pairs):
    if i == j: V[3 * i + i, c] = 1
    else: V[3 * i + j, c] = V[3 * j + i, c] = 1 / np.sqrt(2)
def sym2(U): return V.T @ np.kron(U, U) @ V          # 6x6 rep
# ---- exact twirl on (C^3)^{(x)n} ----
_cache = {}
def perm_data(n):
    if n in _cache: return _cache[n]
    perms = list(itertools.permutations(range(n))); dim = 3 ** n
    idx = np.array(list(itertools.product(range(3), repeat=n)))
    pw = 3 ** np.arange(n - 1, -1, -1)
    maps = np.array([(idx[:, list(p)] @ pw) for p in perms])                 # P_sigma e_j = e_{maps[s][j]}
    def cyc(p):
        seen, c = set(), 0
        for s in range(n):
            if s in seen: continue
            c += 1; t = s
            while t not in seen: seen.add(t); t = p[t]
        return c
    inv = {p: tuple(np.argsort(p)) for p in perms}
    G = np.array([[3.0 ** cyc(tuple(np.array(a)[list(inv[b])])) for b in perms] for a in perms])
    # verify Gram numerically on a few entries
    for a in range(3):
        for b in range(3): assert abs(G[a, b] - np.sum(maps[a] == maps[b])) < 1e-9
    _cache[n] = (perms, maps, np.linalg.pinv(G, rcond=1e-10), dim); return _cache[n]
def twirl(B, n):
    perms, maps, Gp, dim = perm_data(n); ar = np.arange(dim)
    b = np.array([B[m, ar].sum() for m in maps])               # tr(B P_sigma^dagger)
    c = Gp @ b
    T = np.zeros((dim, dim), complex)
    for cs, m in zip(c, maps): T[m, ar] += cs
    return T
def ptrans(A, n, sites):
    t = A.reshape((3,) * (2 * n)); ax = list(range(2 * n))
    for s in sites: ax[s], ax[n + s] = ax[n + s], ax[s]
    return t.transpose(ax).reshape(3 ** n, 3 ** n)
def haar_U3(rng):
    Z = (rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))) / np.sqrt(2); Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / abs(np.diag(R)))
V2, V3 = np.kron(V, V), np.kron(np.kron(V, V), V)
def M2_H(phi):
    p9 = V @ phi; B = np.kron(np.outer(p9, p9.conj()), np.outer(p9, p9.conj()))
    return V2.T @ twirl(B, 4) @ V2
def R_H(phi):
    p9 = V @ phi; w = np.kron(np.kron(p9.conj(), p9), p9)
    A = np.outer(w, w.conj())
    return V3.T @ ptrans(twirl(ptrans(A, 6, [0, 1]), 6), 6, [0, 1]) @ V3
# ---- Hessian group = qutrit Clifford group mod phases ----
Xs = np.roll(np.eye(3), 1, axis=0); Zc = np.diag([1, w3, w3 ** 2])
F = np.array([[w3 ** (i * j) for j in range(3)] for i in range(3)]) / np.sqrt(3); S = np.diag([1, 1, w3])
def canon(U):
    f = U.ravel(); i = np.flatnonzero(abs(f) > 1e-9)[0]; U = U * (abs(f[i]) / f[i])
    return tuple(np.round(np.concatenate([U.real.ravel(), U.imag.ravel()]), 8)), U
grp = {}; k0, I3 = canon(np.eye(3, dtype=complex)); grp[k0] = I3; fr = [I3]
while fr:
    nx = []
    for U in fr:
        for g in (Xs, Zc, F, S):
            k_, W = canon(g @ U)
            if k_ not in grp: grp[k_] = W; nx.append(W)
    fr = nx
HG = list(grp.values()); U6s = [sym2(U) for U in HG]
irr = np.mean([abs(np.trace(U)) ** 2 for U in U6s])       # = 1 iff irreducible on C^6
def kp(v, p):
    o = np.array([1.0 + 0j])
    for _ in range(p): o = np.kron(o, v)
    return o
def orbit(phi): return np.array([U @ phi for U in U6s])
def M2_X(phi):
    Vv = np.array([kp(s, 2) for s in orbit(phi)]); return Vv.T @ Vv.conj() / len(Vv)
def R_X(phi):
    Vv = np.array([np.kron(s.conj(), kp(s, 2)) for s in orbit(phi)]); return Vv.T @ Vv.conj() / len(Vv)
def lmax(A): return float(np.linalg.eigvalsh((A + A.conj().T) / 2)[-1])
def main():
    rng = np.random.default_rng(20261016); out = {"hessian_order": len(HG), "irreducibility_sum_|tr|^2": float(irr)}
    # validate the exact twirl: invariance under random SU(3) and trace preservation
    p0 = rng.normal(size=6) + 1j * rng.normal(size=6); p0 /= np.linalg.norm(p0)
    RH0 = R_H(p0); U = haar_U3(rng); U6 = sym2(U); W3 = np.kron(np.kron(U6.conj(), U6), U6)
    out["twirl_invariance_residual"] = float(np.abs(W3 @ RH0 @ W3.conj().T - RH0).max()); out["trR_H"] = float(np.trace(RH0).real)
    # M_2(H) also equals the partial trace of R_H over the conjugate factor
    out["M2H_vs_ptrace_RH"] = float(np.abs(M2_H(p0) - np.einsum('aiaj->ij', RH0.reshape(6, 36, 6, 36))).max())
    rows = []
    for seed in range(4):
        r_ = np.random.default_rng(100 + seed); x0 = r_.normal(size=12)
        def phi_of(x): p = x[:6] + 1j * x[6:]; return p / np.linalg.norm(p)
        def resid(x):
            Dm = M2_X(phi_of(x)) - M2_H(phi_of(x)); return np.concatenate([Dm.real.ravel(), Dm.imag.ravel()])
        sol = least_squares(resid, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        phi = phi_of(sol.x)
        RX, RH = R_X(phi), R_H(phi)
        M2X, M2H = M2_X(phi), M2_H(phi)
        PX, PH = 6 * lmax(RX), 6 * lmax(RH)
        rows.append(dict(seed=seed, M2_diff=float(np.abs(M2X - M2H).max()), mixed_diff=float(np.abs(RX - RH).max()),
                         P_X=PX, P_H=PH, gap=PX - PH, thmB=6 * lmax(M2X), untuned_gap=6 * lmax(R_X(phi_of(x0))) - 6 * lmax(R_H(phi_of(x0)))))
        print(rows[-1], flush=True)
    out["rows"] = rows
    json.dump(out, open("r8_verify.json", "w"), indent=1, default=float)
    print({k: v for k, v in out.items() if k != "rows"})
if __name__ == "__main__":
    main()
