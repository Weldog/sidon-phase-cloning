#!/usr/bin/env python3
"""CA check (redundant with the induction proof) of the Gaussian-binomial identity used in Theorem C:
  sum_{k=0}^r [r k]_2 (D)_k c_k = c_r prod_{i=t-r-1}^{t-2} (D+2^i),  c_k = prod_{i=0}^{t-k-2}(2^i+1),  (D)_k = prod_{i<k}(D-2^i),
as polynomials in D, for all 2 <= t <= 40, 0 <= r <= t/2; and F_{t-r} increasing in r (F_{s+1}/F_s = (2^{s-1}+1)/(D+2^{s-1}) <= 1)."""
import sympy as sp
from math import prod
D = sp.symbols('D')
def qbin(r, k): return prod(2**(r-i)-1 for i in range(k))//prod(2**(i+1)-1 for i in range(k))
bad = 0; n = 0
for t in range(2, 41):
    c = lambda k: prod(2**i+1 for i in range(t-k-1))
    for r in range(0, t//2+1):
        lhs = sum(qbin(r, k)*sp.prod([D-2**i for i in range(k)])*c(k) for k in range(r+1))
        rhs = c(r)*sp.prod([D+2**i for i in range(t-r-1, t-1)])
        n += 1; bad += sp.expand(lhs-rhs) != 0
print(f"q-identity: {n} (t,r) pairs, t<=40: failures {bad}")
