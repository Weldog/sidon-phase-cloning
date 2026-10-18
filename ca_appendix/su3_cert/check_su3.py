#!/usr/bin/env python3
"""Standalone checker: SU(3)/Sym^2 kill, G = Hessian group (216 elements mod phases), (T,k) = (1,1), d = 6.
Ball arithmetic only (python-flint 0.9.0 / Arb, 256 bits).  numpy is NOT used; no solver, no eigensolver.
Data file su3_data.json: rational fiducial f_q, test vector v_top, mu, congruence matrix Vh, functionals Bas, test
fiducials fa, fb, box directions v1, v2, preconditioner Ainv, box half-width h.  Every float in the data is read as the
exact dyadic rational it denotes; the checks are valid for ANY such data (the data only make them succeed).
Checks:
 (G)  216 words in the generators X,Z,F,S give pairwise projectively distinct unitaries (|tr A^dag B| < 3 certified);
      with |<X,Z,F,S>/scalars| = 216 [cited, qutrit Clifford/Hessian group] the list is the whole group.
 (C)  commutant dimensions: (1/216) sum |chi_{Sym^2 V}(g)|^2 = 4 (ball), SU(3) side (4,0)+(0,2): 2 => dim S = 2.
 (P)  Casimir-Lagrange projectors have traces 60,24,6+6,15,15 on K and 15,6 on Sym^2 V (consistency).
 (R)  R_X exact (ball); R_H = SU(3)-twirl of R_X (Schur, exact); lambda_max(R_X) >= Rayleigh(v_top);
      mu I - R_H > 0 by congruence with Vh + Gershgorin (both Vh^dag Vh and Vh^dag(mu-R_H)Vh diagonally dominant).
 (M)  exact moment matching: Poincare-Miranda for G(s,t) = Ainv * ell(D2(f_q + s v1 + t v2)) on [-h,h]^2
      (D2 homogeneous => face maps are exact quartics; coefficients by Lagrange interpolation), ell injective on S
      (det [ell D2(fa), ell D2(fb)] != 0), gap loss <= 288 * dist / |f_q|.
Prints PASS iff gap_lower - loss > 0 and every check holds."""
import itertools, json, sys
from fractions import Fraction as Fr
from flint import acb, arb, acb_mat, ctx
ctx.prec = 256
D = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'su3_data.json'))
S2 = arb(2).sqrt(); S3 = arb(3).sqrt(); I = acb(0, 1)
def A_(x): return acb(arb(x[0]), arb(x[1]))
def eye(n): return acb_mat([[1 if i == j else 0 for j in range(n)] for i in range(n)])
def dag(A): return A.conjugate().transpose()
def tr(A): return sum((A[i, i] for i in range(A.nrows())), acb(0))
def kron(A, B):
    ra, ca, rb, cb = A.nrows(), A.ncols(), B.nrows(), B.ncols()
    return acb_mat([[A[i // rb, j // cb] * B[i % rb, j % cb] for j in range(ca * cb)] for i in range(ra * rb)])
def occs(d, m):
    out = []
    for c in itertools.combinations_with_replacement(range(d), m):
        n = [0] * d
        for i in c: n[i] += 1
        out.append(n)
    return out
def Wiso(d):   # isometry Sym^2(C^d) -> C^d (x) C^d
    O = occs(d, 2); rows = [[acb(0)] * len(O) for _ in range(d * d)]
    for c, n in enumerate(O):
        idx = [i for i in range(d) for _ in range(n[i])]; perms = set(itertools.permutations(idx)); nr = arb(len(perms)).sqrt()
        for p in perms: rows[p[0] * d + p[1]][c] = acb(1) / nr
    return acb_mat(rows)
W32, W62 = Wiso(3), Wiso(6)
w = (acb(-1) + I * S3) / 2
gens = [acb_mat([[0, 0, 1], [1, 0, 0], [0, 1, 0]]), acb_mat([[1, 0, 0], [0, w, 0], [0, 0, w * w]]),
        acb_mat([[w ** (i * j) / S3 for j in range(3)] for i in range(3)]), acb_mat([[1, 0, 0], [0, 1, 0], [0, 0, w]])]
def projeq_possible(A, B):   # False only if certified |tr(A^dag B)| < 3, i.e. A, B not proportional
    return not (abs(tr(dag(A) * B)) < 3)
grp = [eye(3)]; frontier = [0]
while frontier:            # BFS; a candidate is new unless it may coincide with a listed element (search only)
    new = []
    for a in frontier:
        for g in gens:
            B = g * grp[a]
            if not any(projeq_possible(B, C) for C in grp): grp.append(B); new.append(len(grp) - 1)
    frontier = new
ok_G = len(grp) == 216 and all(abs(tr(dag(grp[i]) * grp[j])) < 3 for i in range(216) for j in range(i))
U6 = [dag(W32) * kron(U, U) * W32 for U in grp]
chi2 = sum((abs((tr(U) ** 2 + tr(U * U)) / 2) ** 2 for U in U6), arb(0)) / 216
ok_C = bool(abs(chi2 - 4) < arb(1) / 4)
gm = []
for a, b in itertools.combinations(range(3), 2):
    E1 = [[acb(0)] * 3 for _ in range(3)]; E1[a][b] = acb(1); E1[b][a] = acb(1); gm.append(acb_mat(E1))
    E2 = [[acb(0)] * 3 for _ in range(3)]; E2[a][b] = -I; E2[b][a] = I; gm.append(acb_mat(E2))
gm += [acb_mat([[1, 0, 0], [0, -1, 0], [0, 0, 0]]), acb_mat([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) * (acb(1) / S3)]
def lift2(g, W, d): return dag(W) * (kron(g, eye(d)) + kron(eye(d), g)) * W
g6 = [lift2(g, W32, 3) for g in gm]; g21 = [lift2(g, W62, 6) for g in g6]
Gk = [kron(-(g.conjugate()), eye(21)) + kron(eye(6), h) for g, h in zip(g6, g21)]
CK = Gk[0] * Gk[0]
for g in Gk[1:]: CK = CK + g * g
CK = CK * (acb(1) / 4); C21 = g21[0] * g21[0]
for g in g21[1:]: C21 = C21 + g * g
C21 = C21 * (acb(1) / 4)
def cas(p, q): return arb(p * p + q * q + p * q + 3 * p + 3 * q) / 3
def lagrange(C, n, lam, others):
    P = eye(n)
    for m in others: P = P * (C - eye(n) * m) * (acb(1) / (lam - m))
    return P
labs = [(4, 2), (3, 1), (2, 0), (0, 4), (1, 2)]
PK = {l: lagrange(CK, 126, cas(*l), [cas(*m) for m in labs if m != l]) for l in labs}
P40 = lagrange(C21, 21, cas(4, 0), [cas(0, 2)]); P02 = lagrange(C21, 21, cas(0, 2), [cas(4, 0)])
QA = PK[(2, 0)] * kron(eye(6), P40); QB = PK[(2, 0)] * kron(eye(6), P02)
want = {(4, 2): 60, (3, 1): 24, (2, 0): 12, (0, 4): 15, (1, 2): 15}
ok_P = all(abs(tr(PK[l]) - want[l]) < arb(10) ** -20 for l in labs) and abs(tr(P40) - 15) < arb(10) ** -20 and abs(tr(P02) - 6) < arb(10) ** -20 \
    and abs(tr(QA) - 6) < arb(10) ** -20 and abs(tr(QB) - 6) < arb(10) ** -20
iota = acb_mat(126, 6)
for j in range(6):
    for i in range(6):
        e = [acb(0)] * 36; e[i * 6 + j] = acb(1); col = dag(W62) * acb_mat([[x] for x in e])
        for c in range(21): iota[i * 21 + c, j] += col[c, 0]
tau0 = QB * iota * dag(iota) * QA; tau = tau0 * (acb(1) / (tr(dag(tau0) * tau0) / 6).real.sqrt())
def twirlK(R):
    out = acb_mat(126, 126)
    for l in labs:
        if l != (2, 0): out = out + PK[l] * (tr(PK[l] * R) / want[l])
    out = out + QA * (tr(QA * R) / 6) + QB * (tr(QB * R) / 6)
    c = tr(dag(tau) * QB * R * QA) / 6
    return out + tau * c + dag(tau) * c.conjugate()
O6 = occs(6, 2)
def sv2(x): return [x[[i for i in range(6) for _ in range(n[i])][0]] * x[[i for i in range(6) for _ in range(n[i])][1]] * (S2 if max(n) == 1 else 1) for n in O6]
def orbit(f, normalise=True):
    if normalise:
        nr = sum((abs(z) ** 2 for z in f), arb(0)).sqrt(); f = [z / nr for z in f]
    fv = acb_mat([[z] for z in f]); A = acb_mat(126, 216); Bm = acb_mat(21, 216)
    for c, U in enumerate(U6):
        x = U * fv; x = [x[i, 0] for i in range(6)]; s = sv2(x)
        for i in range(6):
            for j in range(21): A[i * 21 + j, c] = x[i].conjugate() * s[j]
        for j in range(21): Bm[j, c] = s[j]
    return A, Bm
def D2(f, normalise=True):
    A, Bm = orbit(f, normalise); MX = Bm * dag(Bm) * (acb(1) / 216)
    return MX - P40 * (tr(P40 * MX) / 15) - P02 * (tr(P02 * MX) / 6)
fq = [(Fr(a), Fr(b)) for a, b in D['fiducial']]
F = [acb(arb(a.numerator) / a.denominator, arb(b.numerator) / b.denominator) for a, b in fq]
A, Bm = orbit(F); RX = A * dag(A) * (acb(1) / 216); RH = twirlK(RX)
v = acb_mat([[A_(z)] for z in D['v_top']]); ray = (dag(v) * RX * v)[0, 0].real / (dag(v) * v)[0, 0].real
mu = arb(D['mu']); Vh = acb_mat([[A_(z) for z in row] for row in D['Vh']])
def gersh(Bx): return min((Bx[i, i].real - sum((abs(Bx[i, j]) for j in range(126) if j != i), arb(0)) for i in range(126)), key=lambda x: x.lower())
gB = gersh(dag(Vh) * (eye(126) * mu - RH) * Vh); gV = gersh(dag(Vh) * Vh)
ok_R = bool(gB > 0 and gV > 0); gap = 6 * ray - 6 * mu
Bas = D['Bas']
def ell(M):
    out = []
    for b in Bas:
        s = arb(0)
        for idx in range(441):
            i, j = divmod(idx, 21)
            if b[idx] != 0: s += arb(b[idx]) * M[i, j].real
            if b[441 + idx] != 0: s += arb(b[441 + idx]) * M[i, j].imag
        out.append(s)
    return out
la, lb = ell(D2([A_(z) for z in D['fa']])), ell(D2([A_(z) for z in D['fb']])); det = la[0] * lb[1] - la[1] * lb[0]
ok_inj = bool(det > 0 or det < 0)
v1 = [A_(z) for z in D['v1']]; v2 = [A_(z) for z in D['v2']]; Ai = [[arb(x) for x in r] for r in D['Ainv']]; h = arb(D['h'])
def Gpt(s, t):
    l = ell(D2([F[i] + v1[i] * s + v2[i] * t for i in range(6)], normalise=False))
    return [Ai[0][0] * l[0] + Ai[0][1] * l[1], Ai[1][0] * l[0] + Ai[1][1] * l[1]]
def face(comp, fix_s, val):
    nodes = [arb(j) * h / 2 for j in (-2, -1, 0, 1, 2)]
    ys = [Gpt(val, x)[comp] if fix_s else Gpt(x, val)[comp] for x in nodes]; coef = [arb(0)] * 5
    for i, xi in enumerate(nodes):
        poly = [arb(1)]; den = arb(1)
        for j, xj in enumerate(nodes):
            if j != i:
                poly = [(poly[k - 1] if k > 0 else arb(0)) - xj * (poly[k] if k < len(poly) else arb(0)) for k in range(len(poly) + 1)]; den *= xi - xj
        for k in range(5): coef[k] += ys[i] * poly[k] / den
    tail = sum((abs(coef[p]) * h ** p for p in range(1, 5)), arb(0))
    return (1 if coef[0] - tail > 0 else (-1 if coef[0] + tail < 0 else 0))
sp, sm, tp, tm = face(0, True, h), face(0, True, -h), face(1, False, h), face(1, False, -h)
detA = Ai[0][0] * Ai[1][1] - Ai[0][1] * Ai[1][0]
ok_PM = sp * sm == -1 and tp * tm == -1 and bool(detA > 0 or detA < 0)   # Ainv invertible
nf = sum((abs(z) ** 2 for z in F), arb(0)).sqrt()
dist = h * sum((abs(a) ** 2 + abs(b) ** 2 for a, b in zip(v1, v2)), arb(0)).sqrt() * 2   # |s v1 + t v2| <= h(|v1|+|v2|) <= 2h sqrt(|v1|^2+|v2|^2)
loss = 288 * dist / nf
PASS = ok_G and ok_C and ok_P and ok_R and ok_inj and ok_PM and bool(gap - loss > 0)
print('group_216_distinct', ok_G, '| commutant_dim_G(Sym2V)', chi2.str(8), ok_C, '| projector_traces', ok_P)
print('P_X_lower', (6 * ray).str(15), '| P_H_upper', (6 * mu).str(15), '| gershgorin', gB.str(4), gV.str(4), ok_R)
print('injectivity_det', det.str(6), ok_inj, '| PM face signs s+,s-,t+,t-', sp, sm, tp, tm, ok_PM)
print('gap_lower', gap.str(12), '| matching_loss_upper', loss.str(6), '| PASS', PASS)
