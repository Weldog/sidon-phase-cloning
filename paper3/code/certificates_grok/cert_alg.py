#!/usr/bin/env python3
"""R25 exact, n-uniform certificate for lambda_max(X^{(T,k)}(n)) = F_{t-w}(D), D = 2^n, from the counts in alg_<a>_<b>.json.
Chain (each step labelled in the report):
 (1) X_op = Z^{-1} sum_{T in Sigma(q)} R(T)                                  [GNW moment thm + partial transpose (= HL Lemma 6)]
 (2) R(T)R(T') = D^{d(T,T')} R(T o T'), Sigma closed under o                  [PROVED, elementary; CA-checked]
 (3) spec(X_op) subset spec(L_X), L_X = left mult. on the abstract algebra    [PROVED, all n]
 (4) rho(L_X) = rho(B'(D)) (equitable quotient, Perron-Frobenius)             [PROVED]
 (5) charpoly(B'(D)) = prod_i (lambda - e_i(D)) exactly (closed-form roots), F = some e_i, and F - e_j > 0 for all D >= 2
     (positive-coefficient test after D = 2 + s)  =>  rho(B') = F                [CA-PROVED here, symbolic in D]
 (6) y = P Delta^{-1} u' satisfies L_X y = F y; if phi(y) != 0 then F in spec(X_op)  [PROVED]
     phi(y) != 0 certified by <R(T0), phi(y)>_HS = sum_T y_T D^{dim(T cap T0)} != 0  [CA-PROVED here, symbolic in D]
 => lambda_max(X_op) = F_{t-w}(D) for every n with D=2^n where (5),(6) hold."""
import sys, json, sympy as sp
from math import prod
D, s = sp.symbols('D s', positive=True)
def kappa(m): q, r = divmod(m, 8); return 4*q + [0,0,0,0,1,1,2,3][r]
def F(t): return sp.Integer(prod(2**i+1 for i in range(t-1))) / (D*sp.prod([D+2**i for i in range(t-1)]))
def pos_for_all_D_ge(expr, D0):
    num, den = sp.fraction(sp.factor(sp.together(expr)))
    def sg(p):
        P = sp.Poly(sp.expand(p.subs(D, D0+s)), s); c = set(sp.sign(x) for x in P.coeffs()); return c
    sn, sd = sg(num), sg(den)
    return (sn <= {1} and sd <= {1}) or (sn <= {-1} and sd <= {-1})
def certify(a, b, nmax=40):
    r = json.load(open(f"alg_{a}_{b}.json")); t = r['t']; m = r['orbits']; sizes = r['sizes']
    Zg = D*sp.prod([D+2**i for i in range(t-1)])
    assert sp.expand(sum(c*D**d for d, c in r['Z']) - Zg) == 0
    Bp = sp.zeros(m, m)
    for i, j, d, v in r['C']: Bp[i, j] += sp.Integer(v)*D**d
    Bp = Bp/Zg
    w = r['w']; T_, k_ = b, a-b
    assert w == b + kappa(a-b)
    Fv = F(t-w)
    ns = (Bp - Fv*sp.eye(m)).applyfunc(lambda e: sp.factor(sp.together(e))).nullspace(simplify=True)
    out = dict(a=a, b=b, T=T_, k=k_, t=t, w=w, F=str(sp.factor(Fv)), orbits=m, sizes=sizes, eigspace_dim=len(ns))
    if not ns: out['ok'] = False; return out
    u = ns[0]; u = u / u[[i for i in range(m) if u[i] != 0][0]]
    u = u.applyfunc(lambda e: sp.factor(sp.together(e)))
    assert all(sp.simplify(e) == 0 for e in (Bp*u - Fv*u))
    lam = sp.symbols('lambda')
    ev = Bp.eigenvals(); assert sum(ev.values()) == m
    cp = sp.factor(sp.together(Bp.charpoly(lam).as_expr()))
    prodf = sp.prod([(lam - e)**mu for e, mu in ev.items()])
    assert sp.simplify(sp.together(cp - prodf)) == 0       # exact factorisation of the characteristic polynomial
    dom = []
    for e in ev:
        if sp.simplify(e - Fv) == 0: dom.append(True); continue
        dom.append(pos_for_all_D_ge(Fv - e, 2))
    upos = all(dom) and any(sp.simplify(e - Fv) == 0 for e in ev)   # rho(B'(D)) = F(D) for every real D >= 2
    out['eigs_Bprime'] = [f"{sp.factor(e)} (mult {mu})" for e, mu in ev.items()]
    out['u'] = [str(e) for e in u]; out['rho_eq_F_all_D>=2'] = upos
    # phi(y) != 0 : s_r(D) = sum_j u_j/|O_j| sum_m H[j][r][m] D^m
    sr = []
    for rr in range(m):
        val = sum(u[j]/sizes[j]*sum(c*D**mm for mm, c in r['H'][j][rr]) for j in range(m))
        sr.append(sp.factor(sp.together(val)))
    out['s_r'] = [str(e) for e in sr]
    good_n = []
    for n in range(1, nmax+1):
        Dv = 2**n
        uv = [e.subs(D, Dv) for e in u]; sv = [e.subs(D, Dv) for e in sr]
        good_n.append(n if any(x != 0 for x in sv) else -n)
    # for all n: some s_r is a nonzero rational function; its numerator has finitely many roots -> list powers-of-2 roots
    nz = [e for e in sr if e != 0]
    roots_pow2 = None
    if nz:
        num = sp.fraction(sp.together(nz[0]))[0]
        rts = sp.Poly(num, D).all_roots() if sp.Poly(num, D).degree() > 0 else []
        roots_pow2 = sorted(set(int(x) for x in rts if x.is_Integer and x > 0 and (int(x) & (int(x)-1)) == 0))
    out['phi_y_nonzero_rational_function'] = bool(nz); out['pow2_roots_of_s'] = roots_pow2
    out['certified_n_1_to_%d' % nmax] = good_n
    out['ok'] = upos and all(g > 0 for g in good_n)
    out['all_n'] = upos and bool(nz) and roots_pow2 == []
    out['P_formula'] = f"C(D+{T_}-1,{T_}) * {out['F']}"
    return out
if __name__ == "__main__":
    allres = []
    for ab in sys.argv[1:]:
        a, b = map(int, ab.split(",")); o = certify(a, b); allres.append(o)
        print(f"T={o['T']} k={o['k']} ({o['T']}->{o['T']+o['k']}) t={o['t']} w={o['w']} F_(t-w)={o['F']} orbits={o['orbits']} "
              f"eigdim={o.get('eigspace_dim')} rho(B')=F(all D>=2)={o.get('rho_eq_F_all_D>=2')} phi(y)!=0 fn={o.get('phi_y_nonzero_rational_function')} "
              f"pow2roots={o.get('pow2_roots_of_s')}  ==> CERT n=1..40: {o['ok']}  ALL n: {o.get('all_n')}", flush=True)
        print("    eig(B') =", o.get('eigs_Bprime'), "\n    u' =", o.get('u'), "\n    s_r =", o.get('s_r'))
    tag = "_".join(x.replace(",", "-") for x in sys.argv[1:])
    json.dump(allres, open(f"cert_alg_{tag}.json", "w"), indent=1)
