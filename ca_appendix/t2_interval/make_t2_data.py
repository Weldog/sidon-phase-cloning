# data generator for check_t2.py (float choice of Lagrange multipliers; any lam >= 0 is valid by weak duality)
import json, math
from mpmath import iv
import sys; sys.path.insert(0, '/workspace/r10'); import t2_proof as tp
res = json.load(open('/workspace/r10/t2_proof.json'))
F = math.factorial; data = dict(m0={str(d): m for d, m in res['m0'].items()}, residual=res['S_ge_1_cases'], lam={})
for d, k in res['S_ge_1_cases']:
    N = k + 2; A, B = tp.avgs(d, k); thB = iv.sqrt(2 / B)
    base = float((thB * F(N)).mid) / (N * (N - 1) / 2)
    data['lam']['%d,%d' % (d, k)] = [0.0] + [base * 10 ** (e / 10) for e in range(-30, 21)]
json.dump(data, open('t2_data.json', 'w'), indent=0)
