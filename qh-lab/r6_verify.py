"""Round-6 audit (independent QH-Lab code).
(1) Theorem F: optimal fusion psi_x (x) psi_y -> psi_{x+y} = E_q(K)/d^3.
(2) Grok's 'genuine' SU(2) tests use COHERENT fiducials: refined Theorem B, P <= dim(span psi^{T}) * lmax(M_{T+k}),
    already forces P_X = P_H there (so those passes are corollaries, not evidence).
(3) First genuinely non-trivial non-abelian test of Conjecture A: non-coherent fiducial |1,0> of spin-1
    (= real unit vectors in R^3 under SO(3)); icosahedral orbits vs exact quadrature of the continuous orbit."""
import json, itertools, numpy as np, cvxpy as cp
from math import comb
SOLV = dict(solver="SCS", eps=1e-9, max_iters=400000)
def kp(v, p):
    o = np.array([1.0 + 0j])
    for _ in range(p): o = np.kron(o, v)
    return o
def counterfeit(states, w, T, k):
    S = np.asarray(states, complex); d = S.shape[1]; din, dout = d ** T, d ** (T + k)
    X = sum(wi * np.kron(np.outer(kp(s, T), kp(s, T).conj()).conj(), np.outer(kp(s, T + k), kp(s, T + k).conj())) for s, wi in zip(S, w))
    J = cp.Variable((din * dout, din * dout), hermitian=True)
    p = cp.Problem(cp.Maximize(cp.real(cp.trace(J @ X))), [J >> 0, cp.partial_trace(J, [din, dout], axis=1) == np.eye(din)])
    p.solve(**SOLV); return float(p.value)
def refined_thmB(states, w, T, k):
    S = np.asarray(states, complex)
    VT = np.array([kp(s, T) for s in S]); dimV = np.linalg.matrix_rank(VT, tol=1e-9)
    M = sum(wi * np.outer(kp(s, T + k), kp(s, T + k).conj()) for s, wi in zip(S, w))
    return float(dimV * np.linalg.eigvalsh((M + M.conj().T) / 2)[-1]), int(dimV)
out = {}
# ---------- (1) Theorem F ----------
def fusion(K, q):
    d = len(K); x = np.arange(q); S = np.exp(2j * np.pi * np.outer(x, K) / q) / np.sqrt(d)
    X = np.zeros((d ** 3, d ** 3), complex)
    for a in range(q):
        for b in range(q):
            inp = np.kron(S[a], S[b]); o = S[(a + b) % q]
            X += np.kron(np.outer(inp, inp.conj()).conj(), np.outer(o, o.conj()))
    X /= q * q
    J = cp.Variable((d ** 3, d ** 3), hermitian=True)
    p = cp.Problem(cp.Maximize(cp.real(cp.trace(J @ X))), [J >> 0, cp.partial_trace(J, [d * d, d], axis=1) == np.eye(d * d)])
    p.solve(**SOLV)
    E = sum(1 for a, b, c, e in itertools.product(K, repeat=4) if (a - b - c + e) % q == 0)
    return float(p.value), E / d ** 3
out["thmF"] = [dict(K=K, q=q, sdp=v, formula=f) for K, q in [([0, 1, 3], 7), ([0, 1, 2], 7), ([0, 1, 5, 17], 32), ([0, 1, 2, 3], 32)]
               for v, f in [fusion(K, q)]]
# ---------- SO(3) icosahedral rotation group ----------
phi = (1 + 5 ** 0.5) / 2
def rot(axis, ang):
    a = np.asarray(axis, float); a /= np.linalg.norm(a); K_ = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K_ + (1 - np.cos(ang)) * K_ @ K_
G = [np.eye(3)]; gens = [rot([0, 1, phi], 2 * np.pi / 5), rot([1, 1, 1], 2 * np.pi / 3)]
fr = [np.eye(3)]
while fr:
    nx = []
    for A in fr:
        for g in gens:
            B = g @ A
            if not any(np.allclose(B, C, atol=1e-9) for C in G): G.append(B); nx.append(B)
    fr = nx
assert len(G) == 60
def lines(P):  # dedupe +/- n
    out_ = []
    for p in P:
        if not any(abs(abs(p @ q) - 1) < 1e-9 for q in out_): out_.append(p)
    return np.array(out_)
def leg(l, x): return np.polynomial.legendre.legval(x, [0] * l + [1])
def design_sums(P, ls=(2, 4, 6, 8, 10)):
    C = P @ P.T; return {l: float(np.mean(leg(l, C))) for l in ls}
ico_axis = np.array([0, 1, phi]) / np.linalg.norm([0, 1, phi])
L6 = lines([g @ ico_axis for g in G])                       # 6 five-fold axes (icosahedron, 5-design)
# 9-design orbit: fiducial where the l=6 icosahedral invariant vanishes (scan arc 5-fold -> 3-fold axis, bisect)
e = ico_axis; tri = np.array([1, 1, 1]) / 3 ** 0.5
def n_of(t): v = (1 - t) * e + t * tri; return v / np.linalg.norm(v)
f = lambda t: sum(leg(6, e @ (g @ n_of(t))) for g in G)
ts = np.linspace(0.01, 0.99, 200); vals = [f(t) for t in ts]
i = next(j for j in range(len(ts) - 1) if vals[j] * vals[j + 1] < 0); lo, hi = ts[i], ts[i + 1]
for _ in range(80):
    mid = (lo + hi) / 2
    if f(lo) * f(mid) <= 0: hi = mid
    else: lo = mid
L60 = lines([g @ n_of(lo) for g in G])
out["design_check"] = {"L6": design_sums(L6), "L60": design_sums(L60), "n_L60": len(L60)}
# ---------- (2) coherent spin-1 (|1,1>-type) icosahedral tests: refined Thm B vs Grok P_H ----------
def coherent_spin1(n):  # |n> (x) |n> projected to Sym^2(C^2) basis (|00>,(|01>+|10>)/rt2,|11>)
    th, ph = np.arccos(np.clip(n[2], -1, 1)), np.arctan2(n[1], n[0])
    a, b = np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)
    return np.array([a * a, np.sqrt(2) * a * b, b * b])
coh6 = [coherent_spin1(p) for p in [g @ ico_axis for g in G]]       # 12 directions (both antipodes are distinct coherent states)
coh6 = np.array([c for i_, c in enumerate(coh6) if not any(abs(abs(np.vdot(c, d_)) - 1) < 1e-9 for d_ in coh6[:i_])])
out["coherent_icosa12_(1,1)"] = dict(N=len(coh6), refinedB=refined_thmB(coh6, np.ones(len(coh6)) / len(coh6), 1, 1), P_H=3 / 5)
# ---------- (3) NON-coherent fiducial |1,0>: real unit vectors ----------
def quad(deg):  # exact product rule on S^2 for polynomials of degree <= deg
    nz = deg // 2 + 1; nphi = deg + 1
    z, wz = np.polynomial.legendre.leggauss(nz)
    P, W = [], []
    for zi, wi in zip(z, wz):
        for j in range(nphi):
            p_ = 2 * np.pi * j / nphi; r = np.sqrt(1 - zi * zi)
            P.append([r * np.cos(p_), r * np.sin(p_), zi]); W.append(wi / 2 / nphi)
    return np.array(P), np.array(W)
tests = []
for name, P_, T, k in [("icosa 6 lines (t=5)", L6, 1, 1), ("icosa 6 lines (t=5)", L6, 1, 2),
                       ("icosa 30-line orbit (t=9)", L60, 2, 1), ("icosa 30-line orbit (t=9)", L60, 1, 3)]:
    m, mix = T + k, 2 * T + k
    QP, QW = quad(2 * mix + 2)
    PH = counterfeit(QP.astype(complex), QW, T, k)
    PX = counterfeit(P_.astype(complex), np.ones(len(P_)) / len(P_), T, k)
    MX = sum(np.outer(kp(p.astype(complex), m), kp(p.astype(complex), m)) for p in P_) / len(P_)
    MH = sum(w_ * np.outer(kp(p.astype(complex), m), kp(p.astype(complex), m)) for p, w_ in zip(QP, QW))
    RX = sum(np.outer(kp(p.astype(complex), 2 * T + k), kp(p.astype(complex), 2 * T + k)) for p in P_) / len(P_)
    RH = sum(w_ * np.outer(kp(p.astype(complex), 2 * T + k), kp(p.astype(complex), 2 * T + k)) for p, w_ in zip(QP, QW))
    tests.append(dict(orbit=name, T=T, k=k, moment_match=bool(np.abs(MX - MH).max() < 1e-10),
                      mixed_match=bool(np.abs(RX - RH).max() < 1e-10), P_X=PX, P_H=PH, gap=PX - PH,
                      thmB_H=refined_thmB(QP.astype(complex), QW, T, k)[0]))
    print(tests[-1], flush=True)
out["noncoherent_tests"] = tests
json.dump(out, open("r6_verify.json", "w"), indent=1, default=float)
