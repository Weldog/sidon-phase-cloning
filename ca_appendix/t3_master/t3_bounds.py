#!/usr/bin/env python3
"""R11 Theorem T3 (torus packing inequality at T = 3, all d >= 3, k >= 1): rigorous evaluation of the master bound S3
and of the analytic region bounds B1, B1c, B2, E3.  Exact rationals (fractions) + mpmath interval arithmetic (iv).
No solver.  The proof (qhash_round11.md Sec. 4) shows rho(F) <= S3(d,k) for every multi-member (ii)-family, and
  Region U1  (d >= DU, 1 <= k <= d-3): S3 <= B1(d,k) (k < KS) or B1c(d, KS+3) (k >= KS), both decreasing in d;
  Region U2  (d >= DU, k >= d-2):     S3 <= B2(d), decreasing in d;
  Region F   (3 <= d < DU, N = k+3 >= max(m0(d) d, d+1)): S3 <= E3(d, floor(N/d)), nonincreasing in m;
  Finite     (3 <= d < DU, N < max(m0 d, d+1)): S3(d,k) evaluated here; pairs with S3 >= 1 are listed (they are
             covered by exact clique-LP certificates, ca_appendix/t3_clique_lp).
Usage: python3 t3_bounds.py [DU] [KS]"""
import math, sys, itertools, json
from fractions import Fraction as Fr
from mpmath import iv
iv.dps = 30
F = math.factorial
def Mb(N, m):
    if m <= 0: return 0
    q, s = divmod(N, m); return F(N) // (F(q + 1) ** s * F(q) ** (m - s))
def I(x): return iv.mpf(x.numerator) / x.denominator if isinstance(x, Fr) else iv.mpf(x)
def kappas(d, k):   # exact kappa_X = avg_X / M*,  X = A (3e1), B (2e1+e2), C (e1+e2+e3); g uniform over balanced C
    q, s = divmod(k, d); N = k + 3; Ms = Mb(N, d); out = []
    for bvec in ([3], [2, 1], [1, 1, 1]):
        r = len(bvec); tot = Fr(0)
        for hs in itertools.product((0, 1), repeat=r):
            h = sum(hs)
            if h > s or r - h > d - s: continue
            p = Fr(math.comb(d - r, s - h), math.comb(d, s))
            den = math.prod(F(bv + q + hv) for bv, hv in zip(bvec, hs)) * F(q + 1) ** (s - h) * F(q) ** (d - s - (r - h))
            tot += p * Fr(F(N), den)
        out.append(tot / Ms)
    return out
def s3min(l, N): return math.comb(l, 3) + (l - 1 if N > l else 0) + (1 if N >= l + 2 else 0)
def S3(d, k):
    N = k + 3; Ms = Mb(N, d); kA, kB, kC = kappas(d, k)
    thA, thB, thC = iv.sqrt(1 / I(kA)), iv.sqrt(3 / I(kB)), iv.sqrt(6 / I(kC))
    W = d * iv.sqrt(I(kA)) + d * (d - 1) * iv.sqrt(3 * I(kB)) + math.comb(d, 3) * iv.sqrt(6 * I(kC))
    rA2 = Fr(math.comb(N, 2) * Mb(N - 2, d - 1), Ms); r1 = Fr(N * Mb(N - 1, d - 1), Ms); r0 = Fr(Mb(N, d - 1), Ms)
    VL = thA * I(rA2) + (d - 1) * thB * I(r1 + r0) + math.comb(d - 1, 2) * thC * I(r0)
    Nf = max(Fr(N * Mb(N - 1, l) * (2 * d - l), 2 * l * Ms) for l in range(1, d))
    thmax = max([thA, thB, thC], key=lambda x: x.b)
    VM = thmax * d * I(Nf)
    phB = max([Fr(Mb(N, l), Ms * s3min(l, N)) for l in range(1, min(d - 2, N) + 1)] or [Fr(0)])
    phC = max([Fr(Mb(N, l), Ms * s3min(l, N)) for l in range(1, min(d - 3, N) + 1)] or [Fr(0)])
    VP = d * thA * I(r0) + max([thB * I(phB), thC * I(phC)], key=lambda x: x.b) * math.comb(d + 2, 3)
    return (max([VL, VM], key=lambda x: x.b) + VP) / W
SQ = lambda x: iv.sqrt(I(x))
def U_common(d, N, kA, kB, kC, r0, J):   # W0-version region bound; kappas are LOWER bounds, r0 and 1/J UPPER bounds
    C3 = math.comb(d, 3)
    VL = 1 / (C3 * SQ(6 * kA * kC)) + 2 * (d - 1) / (C3 * SQ(2 * kB * kC)) + I(Fr(3, d) * r0 / kC)
    thmax = max([1 / SQ(kA), SQ(Fr(3) / kB), SQ(Fr(6) / kC)], key=lambda x: x.b)
    return VL, thmax
def B1(d, k, crude=False):   # Region U1 bound (N = k+3 <= d); refined kappas unless crude
    N = k + 3; t = Fr(k, d)
    kA, kB, kC = (Fr(1, 24), Fr(1, 12), Fr(1, 8)) if crude else (max(Fr(1, 6) * (1 - Fr(3, 4) * t), Fr(1, 24)), max(Fr(1, 2) * (1 - Fr(7, 6) * t), Fr(1, 12)), max(1 - Fr(3, 2) * t, Fr(1, 8)))
    VL, thmax = U_common(d, N, kA, kB, kC, Fr(1), None); C3 = math.comb(d, 3)
    VM = thmax * I(Fr(d * d, (N - 1) * C3)) / SQ(6 * kC)
    VP = I(Fr(d, C3)) / SQ(6 * kA * kC) + max([1 / SQ(2 * kB * kC), I(1 / kC)], key=lambda x: x.b) * I(Fr(math.comb(d + 2, 3), C3 * math.comb(N, 3)))
    return max([VL, VM], key=lambda x: x.b) + VP
def B2(d):   # Region U2 bound (N >= d+1)
    kA, kB, kC = Fr(1, 24), Fr(1, 12), Fr(1, 8); C3 = math.comb(d, 3)
    VL = 1 / (C3 * SQ(6 * kA * kC)) + I(Fr(3, 2) * (d - 1)) / (C3 * SQ(2 * kB * kC)) + I(Fr(3, d) * Fr(1, 2) / kC)
    thmax = max([1 / SQ(kA), SQ(Fr(3) / kB), SQ(Fr(6) / kC)], key=lambda x: x.b)
    VM = thmax * I(Fr(d * d, (d - 1) * C3)) / SQ(6 * kC)
    VP = I(Fr(d, 2 * C3)) / SQ(6 * kA * kC) + max([1 / SQ(2 * kB * kC), I(1 / kC)], key=lambda x: x.b) * I(Fr(math.comb(d + 2, 3), C3 * 4 * (math.comb(d - 2, 3) + d - 2)))
    return max([VL, VM], key=lambda x: x.b) + VP
def rho_m(m): return Fr(F(m), m ** m)
def E3(d, m):   # Region F bound for N >= max(m d, d+1); full W with kappa lower bounds (q0 >= m-1)
    kC = Fr(m, m + 1) ** 3; kB = Fr(m ** 3, (m + 1) ** 2 * (m + 2)); kA = Fr(m ** 3, (m + 1) * (m + 2) * (m + 3))
    thA, thB, thC = 1 / SQ(kA), SQ(Fr(3) / kB), SQ(Fr(6) / kC)
    W = d * SQ(kA) + d * (d - 1) * SQ(3 * kB) + math.comb(d, 3) * SQ(6 * kC)
    R = rho_m(m)
    r0 = min(Fr(1, 2), R); r1 = min(Fr(1), R * m); rA2 = min(Fr(1), R * m * m / 2) if m >= 2 else Fr(1)
    VL = thA * I(rA2) + (d - 1) * thB * I(r1 + r0) + math.comb(d - 1, 2) * thC * I(r0)
    u = {d - 1: r1}
    for l in range(1, d - 1): u[l] = min(Fr(2) ** (l + 1 - d), R ** (d - 1 - l))
    Nf = max(Fr(2 * d - l, 2 * l) * u[l] for l in range(1, d))
    thmax = max([thA, thB, thC], key=lambda x: x.b); VM = thmax * d * I(Nf)
    ph = lambda a: max([min(Fr(2) ** (l - d), R ** (d - l)) / (math.comb(l, 3) + l) for l in range(1, d - a + 1)] or [Fr(0)])
    VP = d * thA * I(r0) + max([thB * I(ph(2)), thC * I(ph(3))], key=lambda x: x.b) * math.comb(d + 2, 3)
    return (max([VL, VM], key=lambda x: x.b) + VP) / W
if __name__ == '__main__':
    DU = int(sys.argv[1]) if len(sys.argv) > 1 else 100; KS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    out = {}
    b1 = [float(B1(max(DU, k + 3), k).b) for k in range(1, KS)]; b1c = float(B1(DU, KS, crude=True).b); b2 = float(B2(DU).b)
    out['U1_refined_max'] = max(b1); out['U1_crude'] = b1c; out['U2'] = b2
    print('DU=%d KS=%d  U1 refined max %.4f (k<KS) | U1 crude(N=KS+3) %.4f | U2 %.4f' % (DU, KS, max(b1), b1c, b2), flush=True)
    m0 = {}
    for d in range(3, DU):
        m = 1
        while not E3(d, m).b < 1: m += 1
        m0[d] = m
    out['m0'] = m0; out['E3_max_at_m0'] = max(float(E3(d, m0[d]).b) for d in m0)
    print('m0:', m0, ' max E3(d,m0) %.4f' % out['E3_max_at_m0'], flush=True)
    bad = []; worst = (0, None); n = 0
    for d in range(3, DU):
        for k in range(1, max(m0[d] * d, d + 1) - 3):
            v = S3(d, k); n += 1
            if v.b < 1:
                if v.b > worst[0]: worst = (float(v.b), (d, k))
            else: bad.append([d, k])
    out.update(n_finite=n, S3_ge_1=bad, S3_max_below_1=worst)
    print('finite pairs %d ; S3 >= 1 at %d pairs %s ; max S3 < 1: %.5f at %s' % (n, len(bad), bad, worst[0], worst[1]))
    ok = max(b1) < 1 and b1c < 1 and b2 < 1
    out['ANALYTIC_REGIONS_OK'] = ok; print('analytic regions OK:', ok)
    json.dump(out, open('t3_bounds_DU%d.json' % DU, 'w'), indent=1)
