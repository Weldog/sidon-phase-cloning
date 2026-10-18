#!/usr/bin/env python3
"""Standalone checker for typed clique-LP certificates of the torus packing inequality
      sum_{Delta in F} Lambda(Delta) <= Lambda_0   (all admissible families F with >= 2 members).
Pure Python integers / fractions. No solver, no numpy.  Usage: python3 check_torus_lp.py cert.json
Certificate entry: {"d":d,"T":T,"k":k,"mu":[[sigma_partition, "p/q"], ...]}  (missing sigma => mu = 0).
Checks performed, per case (all exact):
  (1) a_b = sqrt(C(T,b)/Z_b) enclosed by integer square roots at scale 10^30;
      Z_b = sum_{g in C} mult(b+g), computed via elementary symmetric polynomial over the balanced set C;
  (2) every shift type D (|D^+| = n+k, |D^-| = n, 1 <= n <= T) is enumerated independently as a
      nonincreasing integer vector; r_D = Lambda_hi(D)/Lambda0_lo;
  (3) cov_D = sum_{c in B_T, |D-c|_1 <= T+k} mu_{type(c)};  nu_D = max(0, r_D - cov_D);
  (4) bound = sum_c mu_{type(c)} + sum_D nu_D  (orbit-weighted)  must be <= 1;
  (5) n = 0 (singleton) types: unbalanced D >= 0 must have r_D < 1 (redundant with 2610.06723 Thm X.1).
Validity: two members of an admissible family are at l1-distance >= 2(T+k+1), so a ball of radius T+k
about any centre c contains at most one member; hence sum_F Lambda/Lambda0 <= sum_F (cov+nu) <= bound."""
import sys, json, math
from fractions import Fraction as Fr
SC = 10 ** 30
def fact(n): return math.factorial(n)
def multi(v): return fact(sum(v)) // math.prod(fact(x) for x in v)
def compositions(n, d):
    if d == 1: yield (n,); return
    for i in range(n, -1, -1):
        for r in compositions(n - i, d - 1): yield (i,) + r
def nonincr(length, total, hi, lo):            # nonincreasing integer vectors, entries in [lo,hi], given sum
    if length == 0:
        if total == 0: yield ()
        return
    for x in range(min(hi, total - lo * (length - 1)), lo - 1, -1):
        if x * length < total: break
        for r in nonincr(length - 1, total - x, x, lo): yield (x,) + r
def orbit_size(v):
    c = {}
    for x in v: c[x] = c.get(x, 0) + 1
    return fact(len(v)) // math.prod(fact(m) for m in c.values())
def check(case):
    d, T, k = case['d'], case['T'], case['k']; R = T + k; q, s = divmod(k, d)
    mu = {tuple(p): Fr(m) for p, m in case['mu']}
    assert all(x >= 0 for x in mu.values())
    B = list(compositions(T, d))
    def Z(b):
        e = [Fr(1)] + [Fr(0)] * s
        for x in b:
            for j in range(s, 0, -1): e[j] += e[j - 1] / (x + q + 1)
        return Fr(fact(T + k), math.prod(fact(x + q) for x in b)) * e[s]
    cache = {}; alo = []; ahi = []
    for b in B:
        key = tuple(sorted(b))
        if key not in cache:
            a2 = Fr(multi(b)) / Z(b); r = math.isqrt(a2.numerator * SC * SC // a2.denominator); cache[key] = (r, r + 1)
        alo.append(cache[key][0]); ahi.append(cache[key][1])
    def Lam(D, A):
        tot = 0
        for b, a in zip(B, A):
            m = [x + y for x, y in zip(b, D)]
            if min(m) >= 0: tot += a * multi(m)
        return tot
    bal = tuple([q + 1] * s + [q] * (d - s)); L0 = Lam(bal, alo)
    ctype = [tuple(sorted([x for x in c if x], reverse=True)) for c in B]
    bound = sum(mu.get(t, 0) for t in ctype)                       # sum over all centres
    worst_single = Fr(0); nact = 0; nshift = 0
    for n in range(0, T + 1):
        for D in nonincr(d, k, k + n, -n):
            if -sum(x for x in D if x < 0) != n: continue
            r = Fr(Lam(D, ahi), L0); w = orbit_size(D); nshift += w
            if n == 0:
                if D != bal: worst_single = max(worst_single, r)
                continue
            nact += 1
            cov = sum(mu.get(t, 0) for c, t in zip(B, ctype) if sum(abs(x - y) for x, y in zip(D, c)) <= R)
            if r > cov: bound += w * (r - cov)
    ok = bound <= 1 and worst_single < 1
    return bound, worst_single, nact, nshift, ok
if __name__ == '__main__':
    data = json.load(open(sys.argv[1])); allok = True; mx = Fr(0)
    for case in data:
        bd, ws, na, nsh, ok = check(case); allok &= ok; mx = max(mx, bd)
        print('d=%d T=%d k=%d  bound<=%.10f  max_unbal_single<=%.6f  active_types=%d shifts=%d  %s' % (
            case['d'], case['T'], case['k'], math.ceil(bd * 10**10) / 10**10, math.ceil(ws * 10**6) / 10**6, na, nsh, 'PASS' if ok else 'FAIL'), flush=True)
    print('CASES=%d  MAX_BOUND<=%.10f  ALL_PASS=%s' % (len(data), math.ceil(mx * 10**10) / 10**10, allok))
