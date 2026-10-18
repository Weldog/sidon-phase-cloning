"""Independent implementation of Round-5 Lemma R (reduced dual for full phase codebooks over Z_q):
P_opt(T->T+k) = min_z sum_w r_T(w)/z_w  s.t.  d^{-(2T+k)} sum_w r_{T+k}(s+w) z_w <= 1  for all s in Z_q.
Written from the statement only (not Grok's core.py). Exact integer representation counts."""
import itertools, math
import numpy as np, cvxpy as cp

def rep_counts_sparse(K, q, m):
    """{v: r_m(v)} by dict convolution (no length-q arrays; works for astronomically large q)."""
    r = {0: 1}
    for _ in range(m):
        nr = {}
        for v, c in r.items():
            for k_ in K:
                u = (v + int(k_)) % q; nr[u] = nr.get(u, 0) + c
        r = nr
    return r

def rep_counts(K, q, m):
    """Dense r_m (kept for small q / backwards compatibility)."""
    r = np.zeros(q, dtype=np.int64)
    for v, c in rep_counts_sparse(K, q, m).items(): r[v] = c
    return r

def lemma_r(K, q, T, k):
    """Reduced dual, sparse constraint matrix + column scaling. Returns a rigorous upper bound (exactly feasible point
    after rescaling) equal to P_opt up to solver tolerance. Rows: s = v - w for v in supp r_{T+k}, w in supp r_T."""
    import scipy.sparse as sp
    d = len(K); logsig = (2 * T + k) * np.log(d)
    rT, rO = rep_counts_sparse(K, q, T), rep_counts_sparse(K, q, T + k)
    W = sorted(rT); col = {w: i for i, w in enumerate(W)}
    rows, cols, vals, rid = [], [], [], {}
    for v, cv in rO.items():
        for w in W:
            s_ = (v - w) % q
            if s_ not in rid: rid[s_] = len(rid)
            rows.append(rid[s_]); cols.append(col[w]); vals.append(float(cv))
    A = sp.csr_matrix((vals, (rows, cols)), shape=(len(rid), len(W)))
    cs = np.asarray(A.max(axis=0).todense()).ravel(); B = (A @ sp.diags(1.0 / cs)).tocsr()
    obj = np.array([rT[w] for w in W], float) * cs; om = obj.max(); obj = obj / om
    z = cp.Variable(len(W), pos=True)
    prob = cp.Problem(cp.Minimize(cp.sum(cp.multiply(obj, cp.inv_pos(z)))), [B @ z <= 1])
    try: prob.solve(solver="CLARABEL")
    except Exception: prob.solve(solver="SCS", eps=1e-10, max_iters=500000)
    zf = np.maximum(np.asarray(z.value).ravel(), 1e-300); zf = zf / max(1.0, (B @ zf).max())
    return float(np.exp(np.log(np.sum(obj / zf)) + np.log(om) - logsig))

def is_Bh(K, q, h):
    sums = [sum(c) % q for c in itertools.combinations_with_replacement([int(x) for x in K], h)]
    return len(set(sums)) == len(sums)

def torus_K(d, T, k):
    """integer K = b^j, b = 2T+k+1: B_{2T+k} over Z; q chosen so no wrap-around (Python ints, any size)."""
    b = 2 * T + k + 1
    K = [b ** j for j in range(d)]
    return K, 1 << int(2 * (2 * T + k) * K[-1] + 2).bit_length()

def t1_closed_form(d, k):
    """(k+1)!(2d-k)/(2d^{k+1}): torus T=1 value, VALID ONLY FOR k <= d (= N=1 case of arXiv:2610.06723 Eq. 16).
    For k > d it is wrong (even negative); use lemma_r on torus_K instead."""
    if k > d: raise ValueError("closed form valid only for k <= d; see arXiv:2610.06723 Eq. (16)")
    return math.factorial(k + 1) * (2 * d - k) / (2 * d ** (k + 1))

def sidon_kraus_fidelity(states):
    """Theorem-1 primal channel K_j|a> = sqrt((2-delta_aj)/(2d-1)) |e_a+e_j>_Sym2, averaged 1->2 fidelity."""
    S = np.asarray(states); d = S.shape[1]
    def sym(a, j):
        v = np.zeros(d * d, complex)
        if a == j: v[a * d + a] = 1
        else: v[a * d + j] = v[j * d + a] = 1 / np.sqrt(2)
        return v
    Ks = []
    for j in range(d):
        M = np.zeros((d * d, d), complex)
        for a in range(d):
            M[:, a] = np.sqrt((2 - (a == j)) / (2 * d - 1)) * sym(a, j)
        Ks.append(M)
    assert np.allclose(sum(K.conj().T @ K for K in Ks), np.eye(d))
    f = []
    for s in S:
        t = np.kron(s, s)
        f.append(sum(abs(t.conj() @ K @ s) ** 2 for K in Ks))
    return float(np.mean(f))

def gf_log_table(p, h):
    """powers of x in GF(p)[x]/(f) for a primitive monic f of degree h -> {tuple: exponent}."""
    for c in itertools.product(range(p), repeat=h):
        if c[0] == 0: continue
        v = [1] + [0] * (h - 1); table = {tuple(v): 0}
        for n in range(1, p ** h):
            top = v[-1]; v = [0] + v[:-1]
            v = [(v[i] - top * c[i]) % p for i in range(h)]
            if tuple(v) in table: break
            table[tuple(v)] = n
        if len(table) == p ** h - 1: return table
    return None

def bose_chowla(p, h):
    tab = gf_log_table(p, h)
    K = sorted(tab[tuple([c, 1] + [0] * (h - 2))] for c in range(p))
    return np.array(K, dtype=np.int64), p ** h - 1
