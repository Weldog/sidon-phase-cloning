#!/usr/bin/env python3
"""R26 pre-registration (written before any t=8 computation).
General claim (Theorem C): lambda_max(X^{(a,b)}) = S_t(w)/Z_{n,t} = F_{t-w}, w = min(a,b)+kappa(|a-b|), where
  S_t(r) := sum_{k=0}^r [r k]_2 (D)_k c_k,  (D)_k = prod_{i<k}(D-2^i),  c_k = prod_{i=0}^{t-k-2}(2^i+1)
and, for every T' with dim N_L(T') = r, the column sum of L_X equals S_t(r)/Z.
Streaming predictions for t=8: #{T in Sigma : K <= N_R(T)} = c_{dim K} for each K in a chain of isotropic subspaces,
and sum_T D^{dim(N_R(T) cap M_r)} = S_8(r) for r = 0..w."""
import sympy as sp
from math import prod, comb
D = sp.symbols('D')
def kappa(m): q, r = divmod(m, 8); return 4*q + [0,0,0,0,1,1,2,3][r]
def qbin(r, k):
    num = prod(2**(r-i)-1 for i in range(k)); den = prod(2**(i+1)-1 for i in range(k)); return num//den
def c(t, k): return prod(2**i+1 for i in range(t-k-1))
def S(t, r): return sp.expand(sum(qbin(r, k)*sp.prod([D-2**i for i in range(k)])*c(t, k) for k in range(r+1)))
def F(t, D_): return sp.Rational(prod(2**i+1 for i in range(t-1)), D_*prod(D_+2**i for i in range(t-1)))
print("t=8 forms: (a,b) | T->T+k | w | c_k (k=0..w) | S_8(r) polynomials | S_8(w)/Z == F_{8-w}")
for (a, b) in [(7,1),(6,2),(5,3),(8,0)]:
    t = 8; w = min(a,b)+kappa(abs(a-b)); Z = D*sp.prod([D+2**i for i in range(t-1)])
    Fs = sp.Integer(prod(2**i+1 for i in range(t-w-1)))/(D*sp.prod([D+2**i for i in range(t-w-1)]))
    print(f"({a},{b}) | {b}->{a} | w={w} | c={[c(t,k) for k in range(w+1)]} | S={[str(sp.factor(S(t,r))) for r in range(w+1)]} | {sp.simplify(S(t,w)/Z - Fs) == 0}")
print("\nPredicted lambda and P = C(D+T-1,T) lambda (T<=3) at n=3,4,5:")
print("| n | T->T+k | t | w | lambda=F_{t-w} | P |")
print("|---|---|---|---|---|---|")
for n in (3,4,5):
    Dv = 2**n
    for (a, b) in [(7,1),(6,2),(5,3)]:
        T = b; w = min(a,b)+kappa(abs(a-b)); lam = F(8-w, Dv)
        print(f"| {n} | {T}->{a} | 8 | {w} | {lam} | {comb(Dv+T-1,T)*lam} |")
