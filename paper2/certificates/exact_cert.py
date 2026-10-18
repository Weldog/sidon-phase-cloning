#!/usr/bin/env python3
"""[CA] Exact primal+dual certificates for P_E(1,k) (T=1) of finite ensembles with amplitudes in Z[i]/2^(r/2).
Coordinates: monomial basis m_n(psi) = prod psi_i^{n_i} on Sym^{1+k}; s = C^{1/2} m with C = diag(multinomials) (x) I_d,
so X = C^{1/2} M C^{1/2},  M = sum_i p_i (m_out(psi_i) (x) conj(psi_i)) (...)^dagger  (exact Gaussian rationals).
DUAL: Y = lam*I_d; I(x)Y - X >= 0  <=>  lam*C^{-1} - M >= 0 (congruence), checked by exact LDL on the real embedding.
PRIMAL: W = exact kernel basis of lam*C^{-1} - M (so C^{1/2}W spans the lam-eigenspace of X); P = projector onto it;
J = (d/r) P is a Choi operator iff Tr_out P = (r/d) I_d, checked exactly: Tr_out P = sum_o c_o [W G^{-1} W^dag]_{(o,i),(o,i')},
G = W^dag C W.  Then <X,J> = d*lam = tr Y.  So P_E(1,k) = d*lam exactly."""
import sys, json, itertools, math
from fractions import Fraction as Fr
import numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import h2_stabilizer as h2
class G:  # Gaussian rational
    __slots__ = ("a", "b")
    def __init__(s, a, b=0): s.a = Fr(a); s.b = Fr(b)
    def __add__(s, o): return G(s.a + o.a, s.b + o.b)
    def __sub__(s, o): return G(s.a - o.a, s.b - o.b)
    def __mul__(s, o): return G(s.a*o.a - s.b*o.b, s.a*o.b + s.b*o.a)
    def conj(s): return G(s.a, -s.b)
    def inv(s): n = s.a*s.a + s.b*s.b; return G(s.a/n, -s.b/n)
    def iszero(s): return s.a == 0 and s.b == 0
    def __repr__(s): return f"{s.a}+{s.b}i"
ZERO, ONE = G(0), G(1)
def exactify(psi):
    """psi (numpy) -> (u Gaussian-int list, r) with psi = phase * u / 2^(r/2)."""
    i0 = int(np.argmax(np.abs(psi) > 1e-9)); psi = psi*np.exp(-1j*np.angle(psi[i0]))
    m2 = max(abs(x)**2 for x in psi); r = round(-math.log2(m2)); s = psi*2**(r/2)
    u = [complex(round(z.real), round(z.imag)) for z in s]
    assert max(abs(a - b) for a, b in zip(u, s)) < 1e-9 and abs(sum(abs(z)**2 for z in u) - 2**r) < 1e-9
    return [G(int(z.real), int(z.imag)) for z in u], r
def build(states, k):
    d = len(states[0]); occ = h2.occs(d, 1 + k); No = len(occ); N = No*d
    M = [[ZERO]*N for _ in range(N)]; p = Fr(1, len(states))
    for psi in states:
        u, r = exactify(psi)
        mo = []
        for o in occ:
            v = ONE
            for i, e in enumerate(o):
                for _ in range(e): v = v*u[i]
            mo.append(v)
        vec = [mo[a]*u[b].conj() for a in range(No) for b in range(d)]   # out (x) conj(in)
        w = p/Fr(2**(r*(2 + k)))
        nz = [(x, vec[x]) for x in range(N) if not vec[x].iszero()]
        for x, vx in nz:
            for y, vy in nz:
                t = vx*vy.conj(); M[x][y] = M[x][y] + G(t.a*w, t.b*w)
    c = [Fr(math.factorial(1 + k), math.prod(math.factorial(e) for e in o)) for o in occ]
    return M, c, occ, d
def psd_exact(H):
    """Exact PSD test of Hermitian H over Q(i) via real symmetric embedding and LDL without pivoting."""
    n = len(H); R = [[Fr(0)]*(2*n) for _ in range(2*n)]
    for x in range(n):
        for y in range(n):
            a, b = H[x][y].a, H[x][y].b
            R[x][y] = a; R[x+n][y+n] = a; R[x][y+n] = -b; R[x+n][y] = b
    m = 2*n; pos = zero = 0
    for kk in range(m):
        pv = R[kk][kk]
        if pv < 0: return False, pos, zero
        if pv == 0:
            if any(R[kk][j] != 0 for j in range(kk+1, m)): return False, pos, zero
            zero += 1; continue
        pos += 1; rowk = R[kk]
        for i in range(kk+1, m):
            f = R[i][kk]
            if f == 0: continue
            f = f/pv; Ri = R[i]
            for j in range(kk+1, m):
                if rowk[j]: Ri[j] -= f*rowk[j]
    return True, pos, zero
def kernel(H):
    n = len(H); A = [row[:] for row in H]; piv = []; rr = 0
    for col in range(n):
        p = next((i for i in range(rr, n) if not A[i][col].iszero()), None)
        if p is None: continue
        A[rr], A[p] = A[p], A[rr]; iv = A[rr][col].inv(); A[rr] = [x*iv for x in A[rr]]
        for i in range(n):
            if i != rr and not A[i][col].iszero():
                f = A[i][col]; A[i] = [A[i][j] - f*A[rr][j] for j in range(n)]
        piv.append(col); rr += 1
    free = [j for j in range(n) if j not in piv]; basis = []
    for fcol in free:
        v = [ZERO]*n; v[fcol] = ONE
        for i, pc in enumerate(piv): v[pc] = ZERO - A[i][fcol]
        basis.append(v)
    return basis
def matinv(A):
    n = len(A); B = [A[i][:] + [ONE if i == j else ZERO for j in range(n)] for i in range(n)]
    for col in range(n):
        p = next(i for i in range(col, n) if not B[i][col].iszero()); B[col], B[p] = B[p], B[col]
        iv = B[col][col].inv(); B[col] = [x*iv for x in B[col]]
        for i in range(n):
            if i != col and not B[i][col].iszero():
                f = B[i][col]; B[i] = [B[i][j] - f*B[col][j] for j in range(2*n)]
    return [row[n:] for row in B]
def certify(name, states, k, out):
    M, c, occ, d = build(states, k); No = len(occ); N = No*d
    Xf = np.array([[complex(float(M[x][y].a), float(M[x][y].b)) for y in range(N)] for x in range(N)])
    Cs = np.sqrt(np.repeat([float(z) for z in c], d)); lamf = np.linalg.eigvalsh(Cs[:, None]*Xf*Cs[None, :]).max()
    lam = Fr(lamf).limit_denominator(100000)
    H = [[(G(lam/c[x//d]) if x == y else ZERO) - M[x][y] for y in range(N)] for x in range(N)]
    ok, pos, zero = psd_exact(H)
    W = kernel(H); r = len(W)
    Cv = [c[x//d] for x in range(N)]
    Gm = [[sum((W[a][x].conj()*W[b][x]*G(Cv[x]) for x in range(N)), ZERO) for b in range(r)] for a in range(r)]
    Gi = matinv(Gm)
    WGW = lambda x, y: sum((W[a][x]*Gi[a][b]*W[b][y].conj() for a in range(r) for b in range(r) if not W[a][x].iszero() and not W[b][y].iszero()), ZERO)
    TrP = [[sum((WGW(o*d + i, o*d + j)*G(c[o]) for o in range(No)), ZERO) for j in range(d)] for i in range(d)]
    target = Fr(r, d); trok = all((TrP[i][j].a == (target if i == j else 0)) and TrP[i][j].b == 0 for i in range(d) for j in range(d))
    val = d*lam; Dk = math.comb(d + k, 1 + k); univ = Fr(d, Dk)
    res = dict(name=name, T=1, k=k, d=d, lam=str(lam), value=str(val), universal=str(univ), dual_PSD_exact=ok,
               LDL_pos_pivots=pos, LDL_zero_pivots=zero, kernel_dim=r, TrOut_P_equals_r_over_d_I=trok,
               certified=bool(ok and trok and r > 0), value_minus_universal=str(val - univ),
               kernel_basis=[[str(z) for z in v] for v in W])
    print(f"{name} 1->{1+k}: value={val} (= {float(val):.9f}) universal={univ} ({float(univ):.9f}); dual PSD exact={ok} "
          f"(pos {pos}, zero {zero}); kernel dim={r}; Tr_out P=(r/d)I exact={trok}; CERTIFIED={res['certified']}", flush=True)
    out.append(res)
if __name__ == "__main__":
    st1 = h2.stab_states(1); st2 = h2.stab_states(2)
    paulis = {"I": np.eye(2), "X": np.array([[0,1],[1,0]]), "Y": np.array([[0,-1j],[1j,0]]), "Z": np.diag([1,-1])}
    sets = [["ZI","IZ"], ["XI","IX"], ["YI","IY"], ["XZ","ZY"], ["YZ","ZX"]]
    mub = []
    for S in sets:
        ops = [np.kron(paulis[s[0]], paulis[s[1]]) for s in S]
        for v in st2:
            if all(abs(abs(np.vdot(v, o @ v)) - 1) < 1e-9 for o in ops): mub.append(v)
    assert len(mub) == 20
    out = []
    for name, st, ks in [("octahedron (1-qubit stabilizer, 6)", st1, [1, 2, 3]), ("2-qubit MUB (20 stabilizer states)", mub, [1, 2]),
                         ("2-qubit stabilizer (60)", st2, [3])]:
        for k in ks: certify(name, st, k, out)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "certificates.json"), "w"), indent=1)
    print("certificates written to certificates.json")
