#!/usr/bin/env python3
"""R26 streaming exact check of the two counting lemmas behind Theorem C, for any (a,b) (used for t <= 8):
 (L1) #{T in Sigma(q_ab) : K <= N_R(T)} = c_{dim K} = prod_{i=0}^{t-k-2}(2^i+1)   for K in two chains of isotropic subspaces;
 (L2) S(M_r) := sum_T D^{dim(N_R(T) cap M_r)} = S_t(r) = c_r prod_{i=t-r-1}^{t-2}(D+2^i)   (column sum of Z*L_X at level r);
 plus |Sigma| = prod_{i=0}^{t-2}(2^i+1), and an isotropy/stochasticity spot check of every 997th T.
Lagrangians are STREAMED (generator; only the (m-1)-level list is stored), one process, checkpoint JSON every 500k."""
import sys, json, time, random; sys.path.insert(0, '.')
from lagr import Space, rref, span_elems, hyperbolic_basis, lagrangians_hyp, popc, check_T
from math import prod
def kappa(m): q, r = divmod(m, 8); return 4*q + [0,0,0,0,1,1,2,3][r]
def lagr_stream(m):
    """same recursion as lagr.lagrangians_hyp, but the top level is a generator"""
    sub = lagrangians_hyp(m-1)
    def lift(v):
        a, b = v & ((1 << (m-1)) - 1), v >> (m-1); return (a << 1) | (b << 1) << m
    e0, f0 = 1, 1 << m; MK = (1 << m) - 1
    Bf = lambda u, v: (popc((u & MK) & (v >> m)) + popc((v & MK) & (u >> m))) & 1
    Qf = lambda u: popc((u & MK) & (u >> m)) & 1
    Wbasis = [lift(1 << i) for i in range(m-1)] + [lift(1 << (m-1+i)) for i in range(m-1)]
    for Lp in sub:
        L = [lift(v) for v in Lp]
        yield [e0] + L
        duals = []
        for j in range(m-1):
            rows = [(sum(Bf(w, L[i]) << idx for idx, w in enumerate(Wbasis)), 1 if i == j else 0) for i in range(m-1)]
            piv = []
            for r, c in rows:
                for (pr, pc, pb) in piv:
                    if (r >> pb) & 1: r ^= pr; c ^= pc
                if r:
                    pb = r.bit_length()-1; piv = [(pr ^ r, pc ^ c, pbb) if (pr >> pb) & 1 else (pr, pc, pbb) for (pr, pc, pbb) in piv]
                    piv.append((r, c, pb))
            sol = 0
            for (pr, pc, pb) in piv:
                if pc: sol |= 1 << pb
            u = 0
            for idx, w in enumerate(Wbasis):
                if (sol >> idx) & 1: u ^= w
            duals.append(u)
        for phi in range(1 << (m-1)):
            basis = [L[j] ^ (e0 if (phi >> j) & 1 else 0) for j in range(m-1)]
            w = 0
            for j in range(m-1):
                if (phi >> j) & 1: w ^= duals[j]
            basis.append(f0 ^ w ^ (e0 if Qf(w) else 0))
            yield basis
def isotropic_chain(S, want_one, rnd):
    """greedy maximal isotropic basis; if want_one, start with 1_t (requires q(1_t)=0)"""
    t = S.t; basis = []
    if want_one:
        one = (1 << t) - 1; assert S.q(one) == 0; basis = [one]
    while True:
        span = set(span_elems(basis))
        c = [v for v in range(1, 1 << t) if S.q(v) == 0 and v not in span and all(S.B(v, u) == 0 for u in basis)]
        if not c: return basis
        basis.append(rnd.choice(c))
def gf2_rank(vs):
    piv = {}
    r = 0
    for v in vs:
        while v:
            h = v.bit_length()-1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; r += 1; break
    return r
def run(a, b, tag):
    S = Space(a, b); t = S.t; m = t-1; F = S.full; ONE2 = (1 << (2*t)) - 1
    w = min(a, b) + kappa(abs(a-b)); rnd = random.Random(2026)
    pairs = hyperbolic_basis(S); assert len(pairs) == m
    def emb(v):
        u = 0
        for i, (e, f) in enumerate(pairs):
            if (v >> i) & 1: u ^= e
            if (v >> (m+i)) & 1: u ^= f
        return u
    chains = []
    chains.append(isotropic_chain(S, False, rnd))
    if S.q((1 << t)-1) == 0: chains.append(isotropic_chain(S, True, rnd))
    else: chains.append(isotropic_chain(S, False, random.Random(7)))
    for ch in chains: assert len(ch) == w, (ch, w)
    # counts[c][r][d] = #{T : dim(N_R(T) cap M_r) = d}
    counts = [[{} for _ in range(w+1)] for _ in chains]
    n = 0; bad = 0; t0 = time.time(); ck = f"stream_{a}_{b}.json"
    for Lb in lagr_stream(m):
        Tb = [emb(v) for v in Lb] + [ONE2]
        # right defect N_R(T) = {y : (0,y) in T}: eliminate on the x-part (low t bits)
        piv = []; NR = []
        for v in Tb:
            for p in piv:
                if (v & F) and ((v & F) ^ (p & F)) < (v & F): v ^= p
            if v & F: piv.append(v); piv.sort(key=lambda z: -(z & F))
            elif v: NR.append(v >> t)
        rNR = gf2_rank(NR)
        if n % 997 == 0:
            T = rref(Tb)
            if not check_T(S, T): bad += 1
        for ci, ch in enumerate(chains):
            for r in range(w+1):
                d = rNR + r - gf2_rank(NR + ch[:r])
                counts[ci][r][d] = counts[ci][r].get(d, 0) + 1
        n += 1
        if n % 500000 == 0:
            json.dump(dict(a=a, b=b, n=n, partial=True, counts=[[sorted(c.items()) for c in cc] for cc in counts]), open(ck, "w"))
            print(f"  ({a},{b}) {n} streamed [{time.time()-t0:.0f}s]", flush=True)
    # verdicts
    c_ = lambda k: prod(2**i+1 for i in range(t-k-1))
    ok_count = (n == c_(0))
    L1 = all(counts[ci][r].get(r, 0) == c_(r) for ci in range(len(chains)) for r in range(w+1))
    # L2: sum_d counts * D^d == c_r prod_{i=t-r-1}^{t-2}(D+2^i)  -> compare polynomial coefficients exactly
    def polymul(p, q_):
        out = [0]*(len(p)+len(q_)-1)
        for i, x in enumerate(p):
            for j, y in enumerate(q_): out[i+j] += x*y
        return out
    L2 = True
    for ci in range(len(chains)):
        for r in range(w+1):
            lhs = [0]*(r+1)
            for d, cnt in counts[ci][r].items(): lhs[d] += cnt
            rhs = [c_(r)]
            for i in range(t-r-1, t-1): rhs = polymul(rhs, [2**i, 1])
            if lhs != rhs: L2 = False
    res = dict(a=a, b=b, t=t, w=w, n=n, partial=False, sigma_count_ok=ok_count, spotcheck_bad=bad,
               chains=chains, L1_ok=L1, L2_ok=L2, counts=[[sorted(c.items()) for c in cc] for cc in counts], secs=time.time()-t0)
    json.dump(res, open(ck, "w"))
    print(f"{tag} (a,b)=({a},{b}) t={t} w={w}: |Sigma| streamed={n} (expected {c_(0)}) ok={ok_count}; spot-check bad={bad}; "
          f"L1 (#T>=K = c_k) ok={L1}; L2 (S(M_r)=S_t(r)) ok={L2}; chains={chains} [{time.time()-t0:.0f}s]", flush=True)
if __name__ == "__main__":
    for ab in sys.argv[1:]:
        a, b = map(int, ab.split(",")); run(a, b, "STREAM")
