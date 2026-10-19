#!/usr/bin/env python3
"""R25 algebra route.  X = Z^{-1} sum_{T in Sigma(q_ab)} R(T) (HL Lemma 6).  Abstract algebra A: basis e_T, e_T e_T' = D^{d(T,T')} e_{T o T'}.
L_X (left multiplication) is entrywise >= 0; it commutes with left/right multiplication by group elements g (invertible T),
so the G x G orbit partition gives the column-equitable quotient  B'_{ij}(D) = Z^{-1} sum_{T} D^{d(T,T'_j)} [T o T'_j in O_i].
Facts used (see report): rho(L_X) = rho(B') (Perron-Frobenius + equitable quotient);  lambda_max(X_op) <= rho(L_X) for all n,
with equality when R(T) are linearly independent (n >= t-1, GNW Thm 4.3).  Counts c_{ijd} are exact integers: B'(D) for every D.
Also computes the code-state Rayleigh quotient  rho_N(D) = Z^{-1} sum_T D^{dim(T cap (N+N)) - w}  (an exact LOWER bound on lambda_max for every n)."""
import sys, json, random, time, pickle; sys.path.insert(0, '.')
from lagr import *
from multiprocessing import Pool
def kappa(m): q, r = divmod(m, 8); return 4*q + [0,0,0,0,1,1,2,3][r]
G_S = None; G_SG = None; G_REPS = None; G_ORB = None
def _init(a, b, Sg, reps, orb):
    global G_S, G_SG, G_REPS, G_ORB
    G_S = Space(a, b); G_SG = Sg; G_REPS = reps; G_ORB = orb
def _work(chunk):
    out = {}
    for ti in chunk:
        T = G_SG[ti]
        for j, Tj in enumerate(G_REPS):
            R, d = compose(G_S, T, Tj)
            i = G_ORB[R]
            out[(i, j, d)] = out.get((i, j, d), 0) + 1
    return out
def max_isotropic(S):
    best = []
    def rec(basis, start):
        nonlocal best
        if len(basis) > len(best): best = list(basis)
        el = span_elems(basis)
        for v in range(start, 1 << S.t):
            if S.q(v) == 0 and all(S.B(v, u) == 0 for u in basis) and v not in el:
                rec(basis + [v], v+1)
                if len(best) == S.t // 2: return
    rec([], 1); return best
def run(a, b, procs=int(__import__("os").environ.get("PROCS","2")), seed=1):
    S = Space(a, b); t = S.t; t0 = time.time()
    Sg = sigma(S); idx = {T: i for i, T in enumerate(Sg)}
    # closure check (sampled exhaustively for small t)
    # group generators: transpositions inside A and inside B, plus random invertible elements
    gens = []
    def perm_img(p): return [1 << p[i] for i in range(t)]
    for i in range(a-1): p = list(range(t)); p[i], p[i+1] = p[i+1], p[i]; gens.append(graph(S, perm_img(p)))
    if a >= 2: p = list(range(a))[1:] + [0] + list(range(a, t)); gens.append(graph(S, perm_img(p)))
    for i in range(a, t-1): p = list(range(t)); p[i], p[i+1] = p[i+1], p[i]; gens.append(graph(S, perm_img(p)))
    inv = [T for T in Sg if len(left_defect(S, T)) == 1 and len(left_defect(S, tuple(((u & S.full) << t) | (u >> t) for u in T))) == 1]
    rnd = random.Random(seed); gens += rnd.sample(inv, min(4, len(inv)))
    for g in gens: assert g in idx
    # union-find orbits under T -> g T, T -> T g
    par = list(range(len(Sg)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    closure_ok = True
    for i, T in enumerate(Sg):
        for g in gens:
            for R, d in (compose(S, g, T), compose(S, T, g)):
                if R not in idx: closure_ok = False; continue
                assert d == 0
                ri, rj = find(i), find(idx[R])
                if ri != rj: par[ri] = rj
    roots = sorted(set(find(i) for i in range(len(Sg))))
    oid = {r: k for k, r in enumerate(roots)}
    orb = {T: oid[find(i)] for i, T in enumerate(Sg)}
    sizes = [0]*len(roots)
    for T in Sg: sizes[orb[T]] += 1
    reps = [Sg[r] for r in roots]
    print(f"(a,b)=({a},{b}) t={t} |Sigma|={len(Sg)} |inv|={len(inv)} orbits={len(roots)} closure_ok={closure_ok} [{time.time()-t0:.0f}s]", flush=True)
    chunks = [list(range(i, len(Sg), procs*4)) for i in range(procs*4)]
    with Pool(procs, initializer=_init, initargs=(a, b, Sg, reps, orb)) as pool:
        parts = pool.map(_work, chunks)
    C = {}
    for p in parts:
        for kk, v in p.items(): C[kk] = C.get(kk, 0) + v
    # Z polynomial: Z(D) = sum_T D^{dim(T cap diag)}
    Zc = {}
    for T in Sg:
        dd = len([1 for u in span_elems(T) if (u & S.full) == (u >> t)]).bit_length() - 1
        Zc[dd] = Zc.get(dd, 0) + 1
    # code-state Rayleigh: N maximal isotropic for q
    N = max_isotropic(S); w = len(N); NN = set()
    Nel = span_elems(N)
    for x in Nel:
        for y in Nel: NN.add(x | (y << t))
    Rc = {}
    for T in Sg:
        m = sum(1 for u in span_elems(T) if u in NN).bit_length() - 1
        Rc[m - w] = Rc.get(m - w, 0) + 1
    # intersection profile H[j][r] = {m: #{T in O_j : dim(T cap T0_r) = m}} for every orbit rep T0_r  (for the phi(y) != 0 test)
    H = [[{} for _ in reps] for _ in reps]
    for T in Sg:
        j = orb[T]
        for r_, T0 in enumerate(reps):
            m_ = 2*t - len(rref(list(T) + list(T0)))
            H[j][r_][m_] = H[j][r_].get(m_, 0) + 1
    res = dict(a=a, b=b, t=t, H=[[sorted(h.items()) for h in row] for row in H], sigma=len(Sg), orbits=len(roots), sizes=sizes, closure_ok=closure_ok,
               C=[[i, j, d, v] for (i, j, d), v in sorted(C.items())], Z=sorted(Zc.items()), w=w, N=N,
               rayleigh=sorted(Rc.items()), secs=time.time()-t0)
    json.dump(res, open(f"alg_{a}_{b}.json", "w"))
    print(f"   done: w={w} N={N} Z={sorted(Zc.items())} [{time.time()-t0:.0f}s]", flush=True)
    return res
if __name__ == "__main__":
    for ab in sys.argv[1:]:
        a, b = map(int, ab.split(","))
        run(a, b)
