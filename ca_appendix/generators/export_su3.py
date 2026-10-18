# R10 OPTIONAL: rigorous re-certification of the R8 SU(3)/Sym^2 kill with a Gaussian-rational fiducial and ball
# arithmetic (python-flint / Arb, 256-bit balls).  H = SU(3) acting on V = Sym^2(C^3) = (2,0) (d = 6), G = Hessian group
# (216 elements mod phases), (T,k) = (1,1), K = conj(V) (x) Sym^2(V) (dim 126).
#  (1) R_X enclosed exactly (group matrices built from exact generators along BFS words).
#  (2) R_H = Twirl_SU(3)(R_X) EXACTLY by Schur: K = [(4,2)+(3,1)+(2,0)_A] + [(0,4)+(1,2)+(2,0)_B]; isotype projectors are
#      Lagrange polynomials in the Casimir (eigenvalues c(p,q) = (p^2+q^2+pq+3p+3q)/3, all distinct); the two (2,0) copies
#      are linked by the intertwiner tau ~ Q_B N Q_A, N = iota iota^dag, iota(v) = sum_i conj(e_i) (x) (e_i . v)  (U(6)-equivariant).
#  (3) P = 6 lambda_max(R) for X and H (Lemma 1; G irreducible on V).  lambda_max(R_X) >= Rayleigh quotient (ball);
#      lambda_max(R_H) <= mu certified by congruence + Gershgorin: B = V^dag (mu I - R_H) V has all Gershgorin discs in (0,inf)
#      and V^dag V has all discs in (0,inf) (V invertible), so mu I - R_H > 0 (Sylvester).
#  (4) Exact matching: M_2(X) - M_2(H) lies in S = {G-invariant, zero-twirl Hermitian ops on Sym^2 V}, real dim 2
#      (G-commutant dim 4, SU(3)-commutant dim 2).  Poincare-Miranda on the preconditioned map G(s,t) = A^{-1} l(D2(f(s,t))),
#      f(s,t) = f_q + s v1 + t v2, box [-h,h]^2, gives an exactly matched f* in the box; l is injective on S (certified).
#  (5) Gap loss bound 144 * |f* - f_q|.
import itertools, json, sys, time, numpy as np
from flint import acb, arb, acb_mat, ctx
ctx.prec = 256
sys.path.insert(0, '/workspace/r7'); sys.path.insert(0, '/workspace/r8')
def occs(d, m):
    out = []
    for c in itertools.combinations_with_replacement(range(d), m):
        n = np.zeros(d, dtype=int)
        for i in c: n[i] += 1
        out.append(n)
    return np.array(out)
t0 = time.time()
S2 = arb(2).sqrt(); S3 = arb(3).sqrt(); I = acb(0, 1)
def M(rows): return acb_mat(rows)
def kron(A, B):
    ra, ca, rb, cb = A.nrows(), A.ncols(), B.nrows(), B.ncols()
    return acb_mat([[A[i // rb, j // cb] * B[i % rb, j % cb] for j in range(ca * cb)] for i in range(ra * rb)])
def eye(n): return acb_mat([[1 if i == j else 0 for j in range(n)] for i in range(n)])
def dag(A): return A.conjugate().transpose()
def tr(A): return sum((A[i, i] for i in range(A.nrows())), acb(0))
def Wiso(d, m):   # isometry Sym^m(C^d) -> (C^d)^{(x)m}, occs order (m = 2 only)
    O = occs(d, m); rows = [[acb(0)] * len(O) for _ in range(d ** m)]
    for c, n in enumerate(O):
        idx = [i for i in range(d) for _ in range(int(n[i]))]
        perms = set(itertools.permutations(idx)); nrm = arb(len(perms)).sqrt()
        for p in perms: rows[p[0] * d + p[1]][c] = acb(1) / nrm
    return acb_mat(rows)
W32, W62 = Wiso(3, 2), Wiso(6, 2)
# exact generators of the Hessian group (as in r8/su3_tests.hessian_group) and BFS words (float only for identification)
w = (acb(-1) + I * S3) / 2
Xg = M([[0, 0, 1], [1, 0, 0], [0, 1, 0]]); Zg = M([[1, 0, 0], [0, w, 0], [0, 0, w * w]])
Fg = M([[w ** (i * j) / S3 for j in range(3)] for i in range(3)]); Sg = M([[1, 0, 0], [0, 1, 0], [0, 0, w]])
gens = [Xg, Zg, Fg, Sg]
def tonp(A): return np.array([[complex(A[i, j].mid().real), complex(A[i, j].mid().imag)] for i in range(A.nrows()) for j in range(A.ncols())])
def fl(A): return np.array([[complex(float(A[i, j].real.mid()), float(A[i, j].imag.mid())) for j in range(A.ncols())] for i in range(A.nrows())])
def key(U):
    u = U.ravel(); i = np.argmax(np.abs(u) > 1e-9); ph = u[i] / abs(u[i]); return tuple(np.round(u / ph, 6))
gf = [fl(g) for g in gens]; grp = [eye(3)]; grpf = [np.eye(3, dtype=complex)]; seen = {key(grpf[0])}; frontier = [0]
while frontier:
    new = []
    for a in frontier:
        for g, gfl in zip(gens, gf):
            Bf = gfl @ grpf[a]; kb = key(Bf)
            if kb not in seen: seen.add(kb); grp.append(g * grp[a]); grpf.append(Bf); new.append(len(grp) - 1)
    frontier = new
assert len(grp) == 216
U6 = [dag(W32) * kron(U, U) * W32 for U in grp]
# su(3) generators (Gell-Mann) lifted
gm = []
for a, b in itertools.combinations(range(3), 2):
    E1 = [[acb(0)] * 3 for _ in range(3)]; E1[a][b] = acb(1); E1[b][a] = acb(1); gm.append(M(E1))
    E2 = [[acb(0)] * 3 for _ in range(3)]; E2[a][b] = -I; E2[b][a] = I; gm.append(M(E2))
gm += [M([[1, 0, 0], [0, -1, 0], [0, 0, 0]]), M([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) * (acb(1) / S3)]
def lift2(g, W, d): return dag(W) * (kron(g, eye(d)) + kron(eye(d), g)) * W
g6 = [lift2(g, W32, 3) for g in gm]; g21 = [lift2(g, W62, 6) for g in g6]
Gk = [kron(-(g.conjugate()), eye(21)) + kron(eye(6), h) for g, h in zip(g6, g21)]
CK = sum(Gk[1:], Gk[0] * Gk[0] - Gk[0] * Gk[0]) if False else None
CK = Gk[0] * Gk[0]
for g in Gk[1:]: CK = CK + g * g
CK = CK * (acb(1) / 4)
C21 = g21[0] * g21[0]
for g in g21[1:]: C21 = C21 + g * g
C21 = C21 * (acb(1) / 4)
def cas(p, q): return arb(p * p + q * q + p * q + 3 * p + 3 * q) / 3
def lagrange(C, n, lam, others):
    P = eye(n)
    for mu in others: P = P * (C - eye(n) * mu) * (acb(1) / (lam - mu))
    return P
labsK = [(4, 2), (3, 1), (2, 0), (0, 4), (1, 2)]; cK = {l: cas(*l) for l in labsK}
PK = {l: lagrange(CK, 126, cK[l], [cK[m] for m in labsK if m != l]) for l in labsK}
P40 = lagrange(C21, 21, cas(4, 0), [cas(0, 2)]); P02 = lagrange(C21, 21, cas(0, 2), [cas(4, 0)])
QA = PK[(2, 0)] * kron(eye(6), P40); QB = PK[(2, 0)] * kron(eye(6), P02)
iota = acb_mat([[acb(0)] * 6 for _ in range(126)])
for j in range(6):
    for i in range(6):
        e = [acb(0)] * 36; e[i * 6 + j] = acb(1); col = dag(W62) * acb_mat([[x] for x in e])
        for c in range(21): iota[i * 21 + c, j] += col[c, 0]
N = iota * dag(iota); tau0 = QB * N * QA; ctau = (tr(dag(tau0) * tau0) / 6).real.sqrt(); tau = tau0 * (acb(1) / ctau)
dims = {l: tr(PK[l]).real for l in labsK}
def twirlK(R):
    out = acb_mat(126, 126)
    for l in labsK:
        if l == (2, 0): continue
        out = out + PK[l] * (tr(PK[l] * R) / dims[l])
    out = out + QA * (tr(QA * R) / 6) + QB * (tr(QB * R) / 6)
    c = tr(dag(tau) * QB * R * QA) / 6
    return out + tau * c + dag(tau) * c.conjugate()
def sv2(x):   # symvec(x, 2) for x in C^6, occs order
    O = occs(6, 2); out = []
    for n in O:
        idx = [i for i in range(6) for _ in range(int(n[i]))]
        coef = S2 if idx[0] != idx[1] else arb(1); out.append(x[idx[0]] * x[idx[1]] * coef)
    return out
def orbit_data(f, normalise=True):
    if normalise:
        nr = sum((abs(z) ** 2 for z in f), arb(0)).sqrt(); f = [z / nr for z in f]
    fv = acb_mat([[z] for z in f])
    A = acb_mat(126, 216); Bm = acb_mat(21, 216)
    for c, U in enumerate(U6):
        x = U * fv; x = [x[i, 0] for i in range(6)]; s = sv2(x)
        for i in range(6):
            for j in range(21): A[i * 21 + j, c] = x[i].conjugate() * s[j]
        for j in range(21): Bm[j, c] = s[j]
    return A, Bm
def M2pair(f, normalise=True):
    A, Bm = orbit_data(f, normalise); M2X = Bm * dag(Bm) * (acb(1) / 216)
    M2H = P40 * (tr(P40 * M2X) / 15) + P02 * (tr(P02 * M2X) / 6)
    return M2X - M2H
# ---- fiducial: R8 float fiducial -> least-squares polish in float -> Gaussian rationals
fk = json.load(open('/workspace/r8/verify_kill.json'))
f0 = np.array(fk['fiducial_re']) + 1j * np.array(fk['fiducial_im'])
from fractions import Fraction as Fr
DEN = 10 ** 18
fq = [(Fr(round(z.real * DEN), DEN), Fr(round(z.imag * DEN), DEN)) for z in f0]
def ball(fr): return [acb(arb(a.numerator) / a.denominator, arb(b.numerator) / b.denominator) for a, b in fr]
F = ball(fq)
A, Bm = orbit_data(F); RX = A * dag(A) * (acb(1) / 216); RH = twirlK(RX)
RXf, RHf = fl(RX), fl(RH)
lx, Vx = np.linalg.eigh((RXf + RXf.conj().T) / 2); lh, Vh = np.linalg.eigh((RHf + RHf.conj().T) / 2)
print('float P_X', 6 * lx[-1], 'P_H', 6 * lh[-1], 'twirl check |[RH,G]|', max(np.abs(fl(RH * g - g * RH)).max() for g in Gk[:2]), time.time() - t0, flush=True)
v = acb_mat([[acb(float(z.real), float(z.imag))] for z in Vx[:, -1]])
ray = (dag(v) * RX * v)[0, 0].real / (dag(v) * v)[0, 0].real
mu = arb(float(lh[-1]) + 1e-10)
Vb = acb_mat([[acb(float(Vh[i, j].real), float(Vh[i, j].imag)) for j in range(126)] for i in range(126)])
Bc = dag(Vb) * (eye(126) * mu - RH) * Vb; Gv = dag(Vb) * Vb
def gersh_pos(Bx):
    worst = None
    for i in range(126):
        r = sum((abs(Bx[i, j]) for j in range(126) if j != i), arb(0)); lo = Bx[i, i].real - r
        worst = lo if worst is None or lo.lower() < worst.lower() else worst
    return worst
gB, gV = gersh_pos(Bc), gersh_pos(Gv)
cholok = bool(gB > 0 and gV > 0)
gap = 6 * ray - 6 * mu
print('ray', ray.str(15), 'gB', gB.str(5), 'gV', gV.str(5), 'gap', gap.str(12), time.time() - t0, flush=True)
# ---- exact matching by Poincare-Miranda
rng = np.random.default_rng(7)
def D2f(fnp):   # float evaluation
    fnp = fnp / np.linalg.norm(fnp); Xs = [fl(U) @ fnp for U in U6]
    O = occs(6, 2); c = np.array([np.sqrt(2) if (n == 1).sum() == 2 else 1 for n in O])
    Ss = np.array([c * np.prod(np.power(x[None, :], O), axis=1) for x in Xs]); Mx = Ss.T @ Ss.conj() / 216
    P40f, P02f = fl(P40), fl(P02); return Mx - P40f * np.trace(P40f @ Mx).real / 15 - P02f * np.trace(P02f @ Mx).real / 6
U6f = [fl(U) for U in U6]
fl_U = lambda U: U
Ds = [D2f(rng.normal(size=6) + 1j * rng.normal(size=6)) for _ in range(4)]
Vs = np.array([np.concatenate([D.real.ravel(), D.imag.ravel()]) for D in Ds]); u_, s_, wt = np.linalg.svd(Vs); Bas = wt[:2]
def ell_np(D): x = np.concatenate([D.real.ravel(), D.imag.ravel()]); return Bas @ x
def ell_b(D):
    out = []
    for b in Bas:
        s = arb(0)
        for idx in range(441):
            i, j = divmod(idx, 21)
            if b[idx] != 0: s += arb(float(b[idx])) * D[i, j].real
            if b[441 + idx] != 0: s += arb(float(b[441 + idx])) * D[i, j].imag
        out.append(s)
    return out
fa = rng.normal(size=6) + 1j * rng.normal(size=6); fb = rng.normal(size=6) + 1j * rng.normal(size=6)
la = ell_b(M2pair([acb(float(z.real), float(z.imag)) for z in fa])); lb = ell_b(M2pair([acb(float(z.real), float(z.imag)) for z in fb]))
det_inj = la[0] * lb[1] - la[1] * lb[0]
fqn = np.array([float(a) + 1j * float(b) for a, b in fq])
def rnd_dir():
    z = rng.normal(size=6) + 1j * rng.normal(size=6); z -= fqn * np.vdot(fqn, z) / np.vdot(fqn, fqn); return z / np.linalg.norm(z)
v1, v2 = rnd_dir(), rnd_dir()
eps = 1e-6
Jm = np.array([(ell_np(D2f(fqn + eps * v)) - ell_np(D2f(fqn - eps * v))) / (2 * eps) for v in (v1, v2)]).T
Ainv = np.linalg.inv(Jm)
h = 1e-12
# D2 is homogeneous of degree (2,2) in (f, conj f) once normalisation is dropped (the matching set is a cone), so on a
# face the preconditioned map is an exact real polynomial of degree <= 4 in the free parameter; its coefficients are
# recovered exactly (ball arithmetic) by Lagrange interpolation at 5 rational nodes, and the face sign is decided by
# |g_0| > sum_{p>=1} |g_p| h^p.
def Gpt(sval, tval):
    f = [acb(arb(a.numerator) / a.denominator, arb(b.numerator) / b.denominator) + acb(float(p.real), float(p.imag)) * sval + acb(float(q.real), float(q.imag)) * tval
         for (a, b), p, q in zip(fq, v1, v2)]
    l = ell_b(M2pair(f, normalise=False)); return [arb(float(Ainv[0, 0])) * l[0] + arb(float(Ainv[0, 1])) * l[1], arb(float(Ainv[1, 0])) * l[0] + arb(float(Ainv[1, 1])) * l[1]]
def face(comp, fixed_is_s, val):
    nodes = [arb(j) * h / 2 for j in (-2, -1, 0, 1, 2)]
    ys = [Gpt(val, x)[comp] if fixed_is_s else Gpt(x, val)[comp] for x in nodes]
    # Lagrange -> monomial coefficients (exact in balls)
    coef = [arb(0)] * 5
    for i, xi in enumerate(nodes):
        poly = [arb(1)]; den = arb(1)
        for j, xj in enumerate(nodes):
            if j == i: continue
            poly = [(poly[k - 1] if k > 0 else arb(0)) - xj * (poly[k] if k < len(poly) else arb(0)) for k in range(len(poly) + 1)]; den *= (xi - xj)
        for k in range(5): coef[k] += ys[i] * poly[k] / den
    tail = sum((abs(coef[p]) * arb(h) ** p for p in range(1, 5)), arb(0))
    return coef[0], tail, (coef[0] - tail > 0), (coef[0] + tail < 0)
fsp, fsm, ftp, ftm = face(0, True, arb(h)), face(0, True, arb(-h)), face(1, False, arb(h)), face(1, False, arb(-h))
pm = bool(((fsp[2] and fsm[3]) or (fsp[3] and fsm[2])) and ((ftp[2] and ftm[3]) or (ftp[3] and ftm[2])))
g_sp, g_sm, g_tp, g_tm = [fsp[0]], [fsm[0]], [None, ftp[0]], [None, ftm[0]]
print('faces', [(x[0].str(5), x[1].str(3)) for x in (fsp, fsm, ftp, ftm)], flush=True)
dist = h * np.sqrt(2) * 1.000001
loss = 2 * 144 * dist / float(sum(abs(fqn) ** 2) ** 0.5)
cert = dict(group_order=216, fiducial_rational=[[str(a), str(b)] for a, b in fq], P_X_float=float(6 * lx[-1]), P_H_float=float(6 * lh[-1]),
            P_X_lower=ray.str(20) if False else float((6 * ray).lower()), P_H_upper=float((6 * mu).upper()), gershgorin_B_lower=float(gB.lower()),
            gershgorin_VdagV_lower=float(gV.lower()), lambdaH_cert_ok=cholok, certified_gap_lower=float(gap.lower()),
            injectivity_det=[float(det_inj.lower()), float(det_inj.upper())], PM_faces=dict(s_plus=g_sp[0].str(5), s_minus=g_sm[0].str(5), t_plus=g_tp[1].str(5), t_minus=g_tm[1].str(5)),
            poincare_miranda_ok=pm, box_halfwidth=h, gap_loss_bound=loss, sec=round(time.time() - t0, 1))
cert['EXACTLY_MATCHED_KILL'] = bool(cholok and pm and (det_inj > 0 or det_inj < 0) and cert['certified_gap_lower'] - loss > 0)
print(cert)

def cl(z): return [float(z.real), float(z.imag)]
data = dict(fiducial=[[str(a), str(b)] for a, b in fq], v_top=[cl(z) for z in Vx[:, -1]], mu=float(lh[-1]) + 1e-10,
            Vh=[[cl(Vh[i, j]) for j in range(126)] for i in range(126)], Bas=[[float(x) for x in b] for b in Bas],
            fa=[cl(z) for z in fa], fb=[cl(z) for z in fb], v1=[cl(z) for z in v1], v2=[cl(z) for z in v2],
            Ainv=[[float(Ainv[i, j]) for j in range(2)] for i in range(2)], h=h)
json.dump(data, open('su3_data.json', 'w'))

