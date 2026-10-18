#!/usr/bin/env python3
"""R13b coverage audit for Theorem packing (c),(d). Imports the shipped checkers (unmodified) and verifies:
 (1) region bookkeeping: every (d,k) in a large box lies in >= 1 region, the finite sets enumerated by the
     checkers equal the analytic 'finite' definitions, and region boundaries are contiguous (no gap);
 (2) region hypotheses used by the analytic bounds hold on every pair assigned to the region
     (N >= d+1 where needed, q0 >= m-1, k <= d-2 / d-3);
 (3) monotonicity claims, spot-checked numerically far beyond the evaluation points;
 (4) direct evaluation of the master bounds S / S3 on the boundary layers of every analytic region (< 1);
 (5) outward (ceil) decimal rounding of every constant quoted in the paper, from the interval upper endpoints."""
import sys, os, math, json, itertools
# Usage: python3 coverage_check.py [path/to/ca_appendix]   (default: ./ca_appendix)
CA = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else 'ca_appendix')
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(CA, 't2_interval')); sys.path.insert(0, os.path.join(CA, 't3_master'))
import check_t2 as C2, t3_bounds as C3
from mpmath import mpf, ceil as mceil
def up(x, n):   # outward decimal: smallest n-digit decimal >= upper endpoint
    v = x.b if hasattr(x, 'b') else x
    return float(mceil(mpf(v) * 10**n)) / 10**n
out = {}; ok = True
# ---------------- T = 2 ----------------
m0_2 = {3: 4, 4: 3, 5: 3, **{d: 2 for d in range(6, 17)}, 17: 1, 18: 1, 19: 1}
D2 = json.load(open(os.path.join(CA, 't2_interval', 't2_data.json')))
assert all(D2['m0'][str(d)] == m0_2[d] for d in m0_2)
def regions2(d, k):
    R = []
    if d >= 40 and k <= d - 2: R.append('I')
    if d >= 20 and k >= d - 1: R.append('II')
    if 3 <= d <= 19 and k >= (m0_2[d] + 1) * d - 2: R.append("II'")
    if (d <= 39 and k <= d - 2) or (d <= 19 and d - 1 <= k <= (m0_2[d] + 1) * d - 3): R.append('fin')
    return R
fin_code = set([(d, k) for d in range(3, 40) for k in range(1, d - 1)] + [(d, k) for d in range(3, 20) for k in range(d - 1, (m0_2[d] + 1) * d - 2)])
fin_def = {(d, k) for d in range(3, 400) for k in range(1, 1200) if 'fin' in regions2(d, k)}
ok2 = fin_code == fin_def and len(fin_code) == 1021
gaps = [(d, k) for d in range(3, 400) for k in range(1, 1200) if not regions2(d, k)]
ok2 &= not gaps
# hypotheses: II and II' need N=k+2 >= d+1; II' needs m=floor(N/d) >= m0 and m <= q0+1
for d in range(3, 400):
    for k in range(1, 1200):
        R = regions2(d, k); N = k + 2; q0 = k // d; m = N // d
        if 'II' in R or "II'" in R: ok2 &= N >= d + 1
        if "II'" in R: ok2 &= (m >= m0_2[d] + 1) and (m <= q0 + 1)
# monotonicity (numerical sanity beyond the evaluation points)
e1 = [C2.E1(d).b for d in range(40, 2001)]; e2u = [C2.E2u(d).b for d in range(20, 2001)]
mono = all(e1[i + 1] <= e1[i] for i in range(len(e1) - 1)) and all(e2u[i + 1] <= e2u[i] for i in range(len(e2u) - 1))
mono2m = all(C2.E2(d, m + 1).b <= C2.E2(d, m).b for d in range(3, 20) for m in range(1, 60))
ok2 &= mono and mono2m
# boundary layers: direct evaluation of S(d,k) inside analytic regions next to their edges
layer = [(40, k) for k in range(1, 39)] + [(41, k) for k in range(1, 40)] + [(d, d - 1) for d in range(20, 60)] + [(d, d) for d in range(20, 60)] \
      + [(d, (m0_2[d] + 1) * d - 2 + j) for d in range(3, 20) for j in range(0, 3)]
smax = max(C2.S(d, k).b for d, k in layer); ok2 &= smax < 1
out['T2'] = dict(finite_code_equals_definition=fin_code == fin_def, n_finite=len(fin_code), gaps_in_box=len(gaps),
                 box='3<=d<400, 1<=k<1200', monotone_E1_E2u_d_to_2000=mono, monotone_E2_in_m_to_60=mono2m,
                 boundary_layer_pairs=len(layer), boundary_layer_max_S=up(smax, 6),
                 ceil_E1_40=up(C2.E1(40), 6), ceil_E2u_20=up(C2.E2u(20), 6),
                 ceil_max_E2_m0=up(max(C2.E2(d, m0_2[d]).b for d in m0_2), 6),
                 ceil_S_9_6=up(C2.S(9, 6), 6), ceil_res_3_1=up(C2.residual_bound(3, 1, D2['lam']['3,1']), 6))
ok &= ok2; out['T2']['OK'] = ok2
# ---------------- T = 3 ----------------
m0_3 = {3: 4, **{d: 3 for d in range(4, 8)}, **{d: 2 for d in range(8, 15)}, **{d: 1 for d in range(15, 100)}}
B = json.load(open(os.path.join(CA, 't3_master', 't3_bounds_DU100.json')))
ok3 = all(B['m0'][str(d)] == m0_3[d] for d in m0_3)
def regions3(d, k):
    N = k + 3; R = []
    if d >= 100 and k <= d - 3: R.append('U1')
    if d >= 100 and k >= d - 2: R.append('U2')
    if 3 <= d <= 99 and N >= max(m0_3[d] * d, d + 1): R.append('F')
    if 3 <= d <= 99 and N < max(m0_3[d] * d, d + 1): R.append('fin')
    return R
fin3_code = {(d, k) for d in range(3, 100) for k in range(1, max(m0_3[d] * d, d + 1) - 3)}
fin3_def = {(d, k) for d in range(3, 400) for k in range(1, 1200) if 'fin' in regions3(d, k)}
ok3 &= fin3_code == fin3_def and len(fin3_code) == 4774 == B['n_finite']
gaps3 = [(d, k) for d in range(3, 400) for k in range(1, 1200) if not regions3(d, k)]
ok3 &= not gaps3
for d in range(3, 400):
    for k in range(1, 1200):
        R = regions3(d, k); N = k + 3; q0 = k // d; m = N // d
        if 'U2' in R or 'F' in R: ok3 &= N >= d + 1
        if 'U1' in R: ok3 &= N <= d
        if 'F' in R: ok3 &= m >= m0_3[d] and m <= q0 + 1
# monotonicity
b1d = all(C3.B1(d + 1, k).b <= C3.B1(d, k).b for k in range(1, 12) for d in range(max(100, k + 3), 600, 7))
b1c = all(C3.B1(d, k + 1, crude=True).b <= C3.B1(d, k, crude=True).b and C3.B1(d + 1, k, crude=True).b <= C3.B1(d, k, crude=True).b
          for d in range(100, 400, 13) for k in range(12, d - 3, 11))
b2 = [C3.B2(d).b for d in range(100, 1500)]; b2m = all(b2[i + 1] <= b2[i] for i in range(len(b2) - 1))
e3m = all(C3.E3(d, m + 1).b <= C3.E3(d, m).b for d in list(range(3, 30)) + [50, 99] for m in range(1, 40))
ok3 &= b1d and b1c and b2m and e3m
# boundary layers: S3 directly, inside F next to the finite edge, and inside U1/U2 at d = 100..102
lay3 = [(d, max(m0_3[d] * d, d + 1) - 3 + j) for d in range(3, 100) for j in range(0, 2)] + [(100, k) for k in (1, 2, 50, 96, 97, 98, 99, 150)] + [(101, k) for k in (1, 97, 98, 99)]
s3max = max(C3.S3(d, k).b for d, k in lay3); ok3 &= s3max < 1
out['T3'] = dict(finite_code_equals_definition=fin3_code == fin3_def, n_finite=len(fin3_code), gaps_in_box=len(gaps3),
                 box='3<=d<400, 1<=k<1200', monotone_B1_d=b1d, monotone_B1crude_dN=b1c, monotone_B2_d_to_1500=b2m, monotone_E3_m_to_40=e3m,
                 boundary_layer_pairs=len(lay3), boundary_layer_max_S3=up(s3max, 6),
                 ceil_U1_refined_max=up(max(C3.B1(max(100, k + 3), k).b for k in range(1, 12)), 6), ceil_U1_crude=up(C3.B1(100, 12, crude=True), 6),
                 ceil_U2=up(C3.B2(100), 6), ceil_max_E3_m0=up(max(C3.E3(d, m0_3[d]).b for d in m0_3), 6), ceil_S3_7_4=up(C3.S3(7, 4), 6))
ok &= ok3; out['T3']['OK'] = ok3
out['ALL_OK'] = ok
print(json.dumps(out, indent=1)); json.dump(out, open('coverage_check.json', 'w'), indent=1)
