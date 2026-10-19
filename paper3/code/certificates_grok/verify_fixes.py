#!/usr/bin/env python3
"""R26 verification of referee fixes R1-01..R1-03 and of the column-sum lemma against the R25 orbit data.
(F1) R1-01: X^{(b,a)} == S_{a,b} conj(X^{(a,b)}) S_{a,b}^dag (block swap), and the OLD identity X^{(b,a)} == conj X^{(a,b)} fails (n=1,2).
(F2) R1-02: every T in Sigma(q_ab) contains (1_t,1_t) [it is the stochastic vector 1_{2t}], so <1|r(T)|1> = 1.
(F3) R1-03: listed base codes are doubly even with dim kappa(m) (m=4..7), extended Hamming [8,4,4] doubly even dim 4;
     direct-sum + pairing construction is totally isotropic for q_{a,b} with dim min(a,b)+kappa(|a-b|), all a,b with t<=24.
(F4) column sums of B' (R25 alg_*.json): sum_i B'_ij(D) == S_t(level j)/Z for every orbit j (t<=7)."""
import sys, json, numpy as np, sympy as sp
sys.path.insert(0, '.'); sys.path.insert(0, '/workspace/r24')
from lagr import *
from h2_stabilizer import stab_states
from math import prod
from orbit_structure_import import orbit_levels
def kappa(m): q, r = divmod(m, 8); return 4*q + [0,0,0,0,1,1,2,3][r]
# F1
for n in (1, 2):
    st = stab_states(n); D = 2**n
    for (a, b) in ([(2,1),(3,1),(1,2),(3,2)] if n == 1 else [(2,1),(3,1)]):
        def mom(a_, b_):
            E = 0
            for s in st:
                v = np.array([1.0+0j])
                for _ in range(a_): v = np.kron(v, s)
                for _ in range(b_): v = np.kron(v, s.conj())
                E = E + np.outer(v, v.conj())
            return E/len(st)
        Xab, Xba = mom(a, b), mom(b, a)
        da, db = D**a, D**b
        Sw = np.zeros((da*db, da*db))
        for i in range(da):
            for j in range(db): Sw[j*da+i, i*db+j] = 1   # |i>_A|j>_B -> |j>_B|i>_A
        new = np.abs(Xba - Sw @ Xab.conj() @ Sw.T).max(); old = np.abs(Xba - Xab.conj()).max()
        print(f"F1 n={n} (a,b)=({a},{b}): |X^(b,a) - S conj(X^(a,b)) S^T| = {new:.1e}   old identity |X^(b,a) - conj X^(a,b)| = {old:.2e}", flush=True)
# F2
bad = 0; tot = 0
for (a, b) in [(3,1),(2,2),(4,1),(3,2),(5,1),(4,2),(3,3),(6,1),(5,2),(4,3)]:
    S = Space(a, b); t = S.t; one = (1 << (2*t)) - 1
    for T in sigma(S):
        tot += 1
        if rref(list(T) + [one]) != T: bad += 1   # (1_t,1_t) = x=1_t (low bits) | y=1_t (high bits) = 1_{2t}
print(f"F2: (1_t,1_t) in T for all {tot} T (t<=7): violations {bad}", flush=True)
# F3
def de_dim(G):
    v = [int(g, 2) for g in G]; el = span_elems(v); return all(popc(x) % 4 == 0 for x in el), len(rref(v)), min(popc(x) for x in el if x)
base = {4: ["1111"], 5: ["11110"], 6: ["111100", "110011"], 7: ["1110100", "0111010", "0011101"]}
H8 = ["11110000", "00111100", "00001111", "01010101"]
for m, G in base.items():
    de, dim, mw = de_dim(G); print(f"F3 base length {m}: doubly even {de}, dim {dim}, kappa({m}) = {kappa(m)}, min wt {mw}")
de, dim, mw = de_dim(H8); print(f"F3 extended Hamming [8,4]: doubly even {de}, dim {dim}, min wt {mw}")
def max_de_code(m):
    q, r = divmod(m, 8); vecs = []; off = 0
    for _ in range(q): vecs += [int(g, 2) << off for g in H8]; off += 8
    if r >= 4: vecs += [int(g, 2) << off for g in base[r]]
    return vecs
bad3 = 0
for t in range(1, 25):
    for b in range(0, t+1):
        a = t-b; S = Space(a, b); p = min(a, b)
        vecs = [(1 << i) | (1 << (a+i)) for i in range(p)]
        sh = p if a >= b else a + p
        vecs += [v << sh for v in max_de_code(abs(a-b))]
        ok = len(rref(vecs)) == p + kappa(abs(a-b)) and all(S.q(v) == 0 for v in vecs) and all(S.B(u, v) == 0 for u in vecs for v in vecs)
        bad3 += (not ok)
print(f"F3 pairing construction isotropic, dim min(a,b)+kappa(|a-b|), all (a,b) with t<=24: failures {bad3}", flush=True)
# F4
D = sp.symbols('D', positive=True)
def Sx(t, r): return prod(2**i+1 for i in range(t-r-1))*sp.prod([D+2**i for i in range(t-r-1, t-1)])
bad4 = 0; cnt = 0
for (a, b) in [(2,1),(3,1),(4,1),(3,2),(5,1),(4,2),(6,0),(6,1),(5,2),(4,3),(7,0)]:
    r = json.load(open(f"/workspace/r25/alg_{a}_{b}.json")); t = a+b; m = r['orbits']
    lev = orbit_levels(a, b)
    for j in range(m):
        col = sum(sp.Integer(v)*D**d for i, jj, d, v in r['C'] if jj == j); cnt += 1
        if sp.expand(col - Sx(t, lev[j])) != 0: bad4 += 1
    print(f"   F4 (a,b)=({a},{b}) levels {lev}", flush=True)
print(f"F4 column sums of Z*B' == S_t(level) for {cnt} orbit columns (t<=7): failures {bad4}")
