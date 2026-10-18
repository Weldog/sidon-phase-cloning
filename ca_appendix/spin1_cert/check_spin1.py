#!/usr/bin/env python3
"""Standalone checker: SU(2)/SO(3) kill at spin 1 (d = 3), G = octahedral rotation group O (24), (T,k) = (1,1).
Ball arithmetic only (python-flint 0.9.0 / Arb, 256 bits); exact integers for the group checks. No solver.
Cartesian basis: O = signed permutation matrices of det 1, H = SO(3).  K = conj(C^3) (x) Sym^2(C^3) (dim 18).
Checks:
 (G) the 24 matrices are distinct and closed under multiplication (exact).
 (C) (1/24) sum_g chi_{Sym^2}(g)^2 = 3 (exact) => O-commutant on Sym^2 C^3 = V_2+V_0 has dim 3; SO(3)-commutant dim 2
     => S = {O-invariant Hermitian ops with zero SO(3)-twirl} is 1-dimensional (real).
 (Q) R_H by the SO(3) product quadrature (7 x 7 uniform in alpha, gamma; 4-node Gauss-Legendre in cos beta, closed-form
     nodes), exact for all matrix coefficients of spin <= 6 (R_H has spin <= 6: degree 3 in g and 3 in conj g).
 (R) lambda_max(R_X) >= Rayleigh(v_top); mu I - R_H > 0 by ball Cholesky (all pivots certified > 0).
     P = 3 lambda_max(R) for both (Lemma 1; O and SO(3) irreducible on C^3).
 (M) c(f) = Re <B4, M_2(X_f) - M_2(H_f)> has c(f_-) < 0 < c(f_+) (or reversed): a zero f* exists on [f_-, f_+]; since
     c(f_-) != 0 and S is 1-dim, B4 is not orthogonal to S, so c(f*) = 0 => M_2(X_{f*}) = M_2(H_{f*}) exactly.
     Loss: |gap(f*) - gap(f_q)| <= 36 * 2 * max|f_+- - f_q| / |f_q|   (u Lipschitz 3, R Lipschitz 6, P = 3 lambda_max).
PASS iff all checks hold and gap_lower - loss > 0."""
import itertools, json, sys
from fractions import Fraction as Fr
from flint import acb, arb, ctx
ctx.prec = 256
D = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'spin1_data.json'))
def matmul(A, B): return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def det3(M): return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                     + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
G = []
for p in itertools.permutations(range(3)):
    for s in itertools.product((1, -1), repeat=3):
        M = tuple(tuple(s[i] if j == p[i] else 0 for j in range(3)) for i in range(3))
        if det3(M) == 1: G.append(M)
Gs = set(G); ok_G = len(G) == 24 and len(Gs) == 24 and all(matmul(a, b) in Gs for a in G for b in G)
tr = lambda M: sum(M[i][i] for i in range(3))
chis = [Fr(tr(g) ** 2 + tr(matmul(g, g)), 2) for g in G]; ok_C = sum(c * c for c in chis) / 24 == 3
PAIRS = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]; SQ2 = arb(2).sqrt()
def sym2(v): return [v[a] * v[b] * (SQ2 if a != b else 1) for a, b in PAIRS]
def quad():
    r = arb(6) / 5; r = r.sqrt(); xa = (arb(3) / 7 - arb(2) / 7 * r).sqrt(); xb = (arb(3) / 7 + arb(2) / 7 * r).sqrt()
    wa = (18 + arb(30).sqrt()) / 36; wb = (18 - arb(30).sqrt()) / 36; pi = arb.pi(); pts = []
    for x, wt in zip([-xb, -xa, xa, xb], [wb, wa, wa, wb]):
        cb, sb = x, (1 - x * x).sqrt()
        for ia in range(7):
            for ig in range(7):
                a = 2 * pi * ia / 7; g = 2 * pi * ig / 7
                Rz = lambda t: [[t.cos(), -t.sin(), arb(0)], [t.sin(), t.cos(), arb(0)], [arb(0), arb(0), arb(1)]]
                Ry = [[cb, arb(0), sb], [arb(0), arb(1), arb(0)], [-sb, arb(0), cb]]
                pts.append((matmul(matmul(Rz(a), Ry), Rz(g)), wt / 2 / 49))
    return pts
PQ = quad(); ok_Q = abs(sum((w for _, w in PQ), arb(0)) - 1) < arb(10) ** -60
def fid(fr):
    F = [acb(arb(Fr(a).numerator) / Fr(a).denominator, arb(Fr(b).numerator) / Fr(b).denominator) for a, b in fr]
    n = sum((abs(z) ** 2 for z in F), arb(0)).sqrt(); return [z / n for z in F], n
def states(F):
    X = [[sum((g[i][k] * F[k] for k in range(3)), acb(0)) for i in range(3)] for g in G]
    Hs = [[sum((R[i][k] * F[k] for k in range(3)), acb(0)) for i in range(3)] for R, _ in PQ]
    return X, [arb(1) / 24] * 24, Hs, [w for _, w in PQ]
def Rmat(st, ws):
    R = [[acb(0)] * 18 for _ in range(18)]
    for v, w in zip(st, ws):
        s = sym2(v); u = [v[a].conjugate() * s[b] for a in range(3) for b in range(6)]
        for i in range(18):
            for j in range(18): R[i][j] += w * u[i] * u[j].conjugate()
    return R
def M2(st, ws):
    M = [[acb(0)] * 6 for _ in range(6)]
    for v, w in zip(st, ws):
        u = sym2(v)
        for i in range(6):
            for j in range(6): M[i][j] += w * u[i] * u[j].conjugate()
    return M
Fq, nq = fid(D['fiducial']); X, wX, Hs, wH = states(Fq)
RX, RH = Rmat(X, wX), Rmat(Hs, wH)
v = [acb(arb(a), arb(b)) for a, b in D['v_top']]
num = sum((v[i].conjugate() * RX[i][j] * v[j] for i in range(18) for j in range(18)), acb(0)).real
lamX = num / sum((abs(z) ** 2 for z in v), arb(0))
mu = arb(D['mu']); A = [[(mu if i == j else 0) - RH[i][j] for j in range(18)] for i in range(18)]
L = [[acb(0)] * 18 for _ in range(18)]; ok_R = True; minpiv = None
for j in range(18):
    s = (A[j][j] - sum((abs(L[j][k]) ** 2 for k in range(j)), acb(0))).real
    if not s > 0: ok_R = False; break
    minpiv = s if minpiv is None or s.lower() < minpiv.lower() else minpiv; L[j][j] = acb(s.sqrt())
    for i in range(j + 1, 18): L[i][j] = (A[i][j] - sum((L[i][k] * L[j][k].conjugate() for k in range(j)), acb(0))) / L[j][j].real
gap = 3 * lamX - 3 * mu
B4 = [[acb(arb(a), arb(b)) for a, b in row] for row in D['B4']]
def cval(fr):
    F, _ = fid(fr); X, wX, Hs, wH = states(F); MX, MH = M2(X, wX), M2(Hs, wH)
    return sum((B4[i][j].conjugate() * (MX[i][j] - MH[i][j]) for i in range(6) for j in range(6)), acb(0)).real
cm, cp = cval(D['f_minus']), cval(D['f_plus'])
ok_M = bool((cm < 0 and cp > 0) or (cm > 0 and cp < 0))
def q2a(x): x = Fr(x); return arb(x.numerator) / x.denominator
def dist(fr): return sum(((q2a(a) - q2a(a0)) ** 2 + (q2a(b) - q2a(b0)) ** 2 for (a, b), (a0, b0) in zip(fr, D['fiducial'])), arb(0)).sqrt()
delta = max(dist(D['f_minus']), dist(D['f_plus']), key=lambda x: x.upper()); loss = 36 * 2 * delta / nq
MXq, MHq = M2(X, wX), M2(Hs, wH); res = max((abs(MXq[i][j] - MHq[i][j]).upper() for i in range(6) for j in range(6)))
PASS = ok_G and ok_C and ok_Q and ok_R and ok_M and bool(gap - loss > 0)
print('group_ok', ok_G, '| commutant_dim_O(Sym2)=3', ok_C, '| quadrature_weights_sum_1', ok_Q)
print('P_X_lower', (3 * lamX).str(15), '| P_H_upper', (3 * mu).str(15), '| cholesky_ok', ok_R, '| min_pivot', minpiv.str(4))
print('M2_residual_at_f_q <= %.3e' % float(res.mid()), '| c(f-)', cm.str(6), '| c(f+)', cp.str(6), '| sign_change', ok_M)
print('gap_lower', gap.str(12), '| matching_loss_upper', loss.str(6), '| PASS', PASS)
