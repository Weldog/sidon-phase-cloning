#!/usr/bin/env python3
"""Theorem C vs every earlier exact Gram-matrix certificate (R23/R24; independent method: exact kernel + Rump PD on the
stabilizer Gram matrix): lambda_cert == F_{t-w}(2^n), t = 2T+k, w = T + kappa(k)."""
import json
from fractions import Fraction as Fr
from math import prod
def kappa(m): q, r = divmod(m, 8); return 4*q + [0,0,0,0,1,1,2,3][r]
def F(t, D): return Fr(prod(2**i+1 for i in range(t-1)), D*prod(D+2**i for i in range(t-1)))
rows = []
for f in ["/workspace/r23/certs_lambda.json", "/workspace/r23/certs_extra.json", "/workspace/r24/certs_r24.json", "/workspace/r24/certs_modular.json"]:
    d = json.load(open(f)); d = d if isinstance(d, list) else [d]
    for r in d:
        if r.get('certified'): rows.append((r['n'], r['T'], r['k'], Fr(r['lam']), f.split('/')[-2]+'/'+f.split('/')[-1]))
bad = 0
print("| n | T->T+k | t | w | certified lambda (source) | F_{t-w} | equal |"); print("|---|---|---|---|---|---|---|")
for n, T, k, lam, src in sorted(set(rows)):
    t = 2*T+k; w = T+kappa(k); f = F(t-w, 2**n); bad += lam != f
    print(f"| {n} | {T}->{T+k} | {t} | {w} | {lam} ({src}) | {f} | {lam == f} |")
print(f"\nmismatches: {bad} of {len(set(rows))}")
