# R11 generator (uses HiGHS LP -- NOT part of the CA checker). Produces typed clique-LP certificates
# mu_sigma (exact rationals) for the torus packing inequality, one JSON per case set.
import sys, json, math, itertools, time
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog

def partitions(n, maxlen=None, maxpart=None):
    if maxpart is None: maxpart = n
    if maxlen is None: maxlen = n
    if n == 0: yield (); return
    if maxlen == 0: return
    for p in range(min(n, maxpart), 0, -1):
        for r in partitions(n - p, maxlen - 1, p): yield (p,) + r
def comps(n, d):
    if d == 1: yield (n,); return
    for i in range(n + 1):
        for r in comps(n - i, d - 1): yield (i,) + r
def mult(v): return math.factorial(sum(v)) // math.prod(math.factorial(x) for x in v)
def orbit(d, *parts):
    used = sum(len(p) for p in parts); den = math.factorial(d - used)
    for p in parts:
        for c in set(p): den *= math.factorial(p.count(c))
    return math.factorial(d) // den
SC = 10 ** 40
def Zb(b, d, k):
    q, s = divmod(k, d); N = sum(b) + k
    base = Fr(math.factorial(N), math.prod(math.factorial(x + q) for x in b))
    e = [Fr(1)] + [Fr(0)] * s
    for x in b:
        w = Fr(1, x + q + 1)
        for j in range(s, 0, -1): e[j] += e[j - 1] * w
    return base * e[s]

def gen(d, T, k):
    t0 = time.time(); R = T + k; q, s = divmod(k, d); bal = tuple([q + 1] * s + [q] * (d - s))
    Bs = list(comps(T, d)); zc = {}
    ahi = {}; alo = {}
    for b in Bs:
        key = tuple(sorted(b))
        if key not in zc:
            a2 = Fr(mult(b)) / Zb(b, d, k); r = math.isqrt(a2.numerator * SC * SC // a2.denominator)
            zc[key] = (r, r + 1)
        alo[b], ahi[b] = zc[key]
    def lam(D, A):
        t = 0
        for b in Bs:
            m = tuple(x + y for x, y in zip(b, D))
            if min(m) >= 0: t += A[b] * mult(m)
        return t
    L0lo = lam(bal, alo)
    sigmas = list(partitions(T, d)); si = {p: i for i, p in enumerate(sigmas)}
    act = []
    for n in range(1, T + 1):
        for ng in partitions(n, d - 1):
            for ps in partitions(n + k, d - len(ng)):
                if len(ng) + len(ps) > d: continue
                D = tuple(list(ps) + [-x for x in ng] + [0] * (d - len(ps) - len(ng)))
                M = [0] * len(sigmas)
                for c in Bs:
                    if sum(abs(x - y) for x, y in zip(D, c)) <= R: M[si[tuple(sorted([x for x in c if x], reverse=True))]] += 1
                act.append((orbit(d, ps, ng), M, Fr(lam(D, ahi), L0lo)))
    ns, nt = len(sigmas), len(act)
    cvec = np.array([orbit(d, p) for p in sigmas] + [a[0] for a in act], float)
    Aub = np.zeros((nt, ns + nt)); bub = np.zeros(nt)
    for i, (sz, M, r) in enumerate(act):
        Aub[i, :ns] = -np.array(M, float); Aub[i, ns + i] = -1; bub[i] = -float(r)
    sc = cvec.max(); res = linprog(cvec / sc, A_ub=Aub, b_ub=bub, bounds=(0, None), method='highs')
    mu = [Fr(max(0.0, x)).limit_denominator(10 ** 12) * Fr(1000000001, 1000000000) for x in res.x[:ns]]
    return dict(d=d, T=T, k=k, lp_float=float(res.fun * sc),
                mu=[[list(p), str(m)] for p, m in zip(sigmas, mu) if m], sec=round(time.time() - t0, 2))

if __name__ == '__main__':
    name = sys.argv[1]
    if name == 'clique_lp_d16_Tk5': cases = [(d, T, k) for d in range(2, 17) for T in range(1, 5) for k in range(1, 6 - T)]
    elif name == 'd3_T60': cases = [(3, T, k) for T in range(1, 41) for k in (1, 2, 3)] + [(3, T, 1) for T in range(41, 61)]
    elif name[0] == 'T':   # T<t>_dlo_dhi_kmax
        t, dlo, dhi, km = name.split('_'); cases = [(d, int(t[1:]), k) for d in range(int(dlo), int(dhi) + 1) for k in range(1, int(km) + 1)]
    out = []; fn = 'cert_%s.json' % name
    for c in cases:
        r = gen(*c); out.append(r); print(c, round(r['lp_float'], 6), r['sec'], flush=True)
        json.dump(out, open(fn, 'w'))
