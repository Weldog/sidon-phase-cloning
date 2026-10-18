#!/usr/bin/env python3
"""Standalone checker for the computational steps of Theorem T2 (torus packing inequality, T = 2, all d >= 3, k >= 1).
Exact rationals (fractions) + mpmath interval arithmetic (iv, 30 digits). No solver.  Usage: python3 check_t2.py t2_data.json
Steps (numbers refer to the proof in qhash_round11.md, Sec. 2.1):
 [C1] E1(40) < 1                       (Region I,  k <= d-2, d >= 40; E1 decreasing in d -- analytic)
 [C2] E2u(20) < 1                      (Region II, k >= d-1, d >= 20; E2u decreasing in d -- analytic)
 [C3] m0(d) minimal with E2(d,m0) < 1, 3 <= d <= 19 (Region II'; E2 nonincreasing in m -- analytic)
 [C4] S(d,k) < 1 on the finite remainder (superset used: d<=39, k<=d-2  and  d<=19, d-1<=k<=(m0+1)d-3), except
      the listed residual pairs
 [C5] residual pairs: Lagrangian bound max(rhoA, rhoB) < 1 with the multipliers lam in the data file.
avg_A, avg_B are computed as exact rationals from the closed form and cross-checked by brute-force enumeration of the
balanced set C for d <= 7."""
import json, math, sys, itertools
from fractions import Fraction as Fr
from mpmath import iv
iv.dps = 30
F = math.factorial
def Mb(N, m):
    if m <= 0: return 0
    q, s = divmod(N, m); return F(N) // (F(q + 1) ** s * F(q) ** (m - s))
def mult(v): return F(sum(v)) // math.prod(F(x) for x in v)
def avgs_exact(d, k):
    q, s = divmod(k, d); N = k + 2
    def m(touched, nh, nl): return Fr(F(N), math.prod(F(x) for x in touched) * F(q + 1) ** nh * F(q) ** nl)
    A = Fr(s, d) * m([q + 3], s - 1, d - s) if s >= 1 else Fr(0)
    A += Fr(d - s, d) * m([q + 2], s, d - s - 1)
    P = d * (d - 1); B = Fr(0)
    if s >= 2: B += Fr(s * (s - 1), P) * m([q + 2, q + 2], s - 2, d - s)
    if s >= 1 and d - s >= 1: B += Fr(2 * s * (d - s), P) * m([q + 2, q + 1], s - 1, d - s - 1)
    if d - s >= 2: B += Fr((d - s) * (d - s - 1), P) * m([q + 1, q + 1], s, d - s - 2)
    return A, B
def avgs_brute(d, k):
    q, s = divmod(k, d); C = set(itertools.permutations([q + 1] * s + [q] * (d - s)))
    bA = [2] + [0] * (d - 1); bB = [1, 1] + [0] * (d - 2)
    return (Fr(sum(mult([x + y for x, y in zip(bA, g)]) for g in C), len(C)), Fr(sum(mult([x + y for x, y in zip(bB, g)]) for g in C), len(C)))
def I(x): return iv.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else iv.mpf(x)
def consts(d, k):
    A, B = avgs_exact(d, k); A, B = I(A), I(B)
    return 1 / iv.sqrt(A), iv.sqrt(2 / B), d * iv.sqrt(A) + iv.mpf(d * (d - 1) // 2) * iv.sqrt(2 * B)
def S(d, k):
    N = k + 2; thA, thB, W = consts(d, k)
    phi = max([Fr(Mb(N, l), l * (l - 1) // 2 + (1 if N > l else 0)) for l in range(1, min(d - 2, N) + 1)] or [Fr(0)])
    num = thA * N * Mb(N - 1, d - 1) + (d - 1) * thB * Mb(N, d - 1) + thB * I(phi) * (d * (d + 1) // 2) + d * thA * Mb(N, d - 1)
    return num / W
def r0(m): return Fr(F(m), m ** m)
def E1(d):
    t3 = max([Fr(1, 3) / (1 - Fr(1, d)), Fr(1, 6) / (1 - Fr(2, d)), Fr(2, 10), Fr(8, (d // 2 + 1) * (d // 2 + 2))])
    return I(Fr(8, d)) + I(Fr(d + 1, d - 1) * t3) + I(Fr(2, d - 1) * (1 + Fr(1, d))) * iv.sqrt(12)
def E2u(d): return I(Fr(8, d) + Fr(d + 1, d - 1) / ((d - 2) * (d - 3) // 2 + 1)) + I(Fr(2, d - 1) * (1 + Fr(1, d))) * iv.sqrt(12)
def E2(d, m):
    kB = Fr(m, m + 1) ** 2; kA = Fr(m * m, (m + 1) * (m + 2))
    rr0 = min(Fr(1), r0(m)); rr1 = min(Fr(1), r0(m) * m); rphi = min(Fr(1, 4 * ((d - 2) * (d - 3) // 2 + 1)), r0(m) ** 2)
    c = iv.sqrt(I(2 * kA * kB))
    return I(Fr(2, d) * rr0 / kB + Fr(d + 1, d - 1) * rphi / kB) + (I(Fr(2, d - 1) * rr0) + I(Fr(2, d * (d - 1)) * rr1)) / c
def parts(n, maxlen, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    if maxlen == 0: return
    for p in range(min(n, maxpart), 0, -1):
        for r in parts(n - p, maxlen - 1, p): yield (p,) + r
def s2(sh): return len(sh) * (len(sh) - 1) // 2 + sum(1 for p in sh if p >= 2)
def bestP(th, avail, N, lam):            # max(0, max_P theta*mult(P) - lam*s2(P)) over shapes with <= avail parts
    b = iv.mpf(0)
    if avail <= 0: return b
    for sh in parts(N, avail):
        v = th * mult(sh) - lam * s2(sh)
        if v.b > b.b: b = v
    return b
def lagr(items, cap, N, lams):           # weak duality: valid for every lam >= 0; take the best upper end
    return min((lam * cap + sum((c * bestP(th, av, N, lam) for c, th, av in items if c > 0), iv.mpf(0)) for lam in lams), key=lambda x: x.b)
def residual_bound(d, k, lams):
    N = k + 2; thA, thB, W = consts(d, k); lams = [iv.mpf(x) for x in lams]; assert all(l.a >= 0 for l in lams)
    rA = lagr([(d, thA, d - 1), (d * (d - 1) // 2, thB, d - 2)], d * (d + 1) // 2, N, lams) / W
    rB = iv.mpf(0)
    for sh in parts(k + 1, d - 1):       # light member Q - e_i, shape of Q
        t = len(sh); o = d - 1 - t
        V1 = mult(sh) * N * (thA + thB * (o + sum((iv.mpf(1) / (p + 1) for p in sh), iv.mpf(0))))
        items = [(t, thA, d - t), (o, thA, d - t - 1), (t * (t - 1) // 2, thB, d - t), (t * o, thB, d - t - 1), (o * (o - 1) // 2, thB, d - t - 2)]
        v = (V1 + lagr(items, (d - t) * (d - t + 1) // 2, N, lams)) / W
        if v.b > rB.b: rB = v
    return max(rA.b, rB.b)
if __name__ == '__main__':
    D = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 't2_data.json')); ok = True
    for d in range(3, 8):
        for k in range(1, 3 * d):
            assert avgs_exact(d, k) == avgs_brute(d, k)
    print('[avg] closed form == brute force for d<=7, k<3d: True')
    e1, e2u = E1(40), E2u(20); ok &= e1.b < 1 and e2u.b < 1
    print('[C1] E1(40) <= %.6f   [C2] E2u(20) <= %.6f' % (float(e1.b), float(e2u.b)))
    print('[C1/C2 sanity] max E1(d), d=40..400: %.6f ; max E2u(d), d=20..400: %.6f' % (max(float(E1(d).b) for d in range(40, 401)), max(float(E2u(d).b) for d in range(20, 401))))
    m0 = {}
    for d in range(3, 20):
        m = 1
        while not E2(d, m).b < 1: m += 1
        m0[d] = m
    ok &= all(m0[d] == D['m0'][str(d)] for d in m0)
    print('[C3] m0 =', m0, ' max E2(d,m0) <= %.6f' % max(float(E2(d, m0[d]).b) for d in m0))
    fin = sorted(set([(d, k) for d in range(3, 40) for k in range(1, d - 1)] + [(d, k) for d in range(3, 20) for k in range(d - 1, (m0[d] + 1) * d - 2)]))
    worst = (0, None); need = []
    for d, k in fin:
        v = S(d, k)
        if v.b < 1:
            if v.b > worst[0]: worst = (float(v.b), (d, k))
        else: need.append([d, k])
    ok &= need == D['residual']
    print('[C4] finite pairs: %d ; S<1 except %d residual pairs %s ; max S below 1: %.6f at %s' % (len(fin), len(need), need, worst[0], worst[1]))
    mx = 0
    for d, k in need:
        b = residual_bound(d, k, D['lam']['%d,%d' % (d, k)]); mx = max(mx, float(b)); ok &= b < 1
        print('[C5] (d,k)=(%d,%d)  rho <= %.6f' % (d, k, float(b)))
    print('[C5] max residual bound %.6f' % mx)
    print('THEOREM T2 COMPUTATIONAL STEPS: %s ; overall sup bound <= %.4f' % ('PASS' if ok else 'FAIL', max(worst[0], mx, float(e1.b), float(e2u.b), max(float(E2(d, m0[d]).b) for d in m0))))
