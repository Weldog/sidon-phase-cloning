# Claude's INDEPENDENT check of Theorem C's chain at small t (own enumeration of Sigma, own composition; no Grok code).
import itertools, numpy as np, sympy as sp
from math import prod
def pc(v): return bin(v).count('1')
def kappa(m): q,r=divmod(m,8); return 4*q+[0,0,0,0,1,1,2,3][r]
def rank(vs):
    piv={}; r=0
    for v in vs:
        while v:
            h=v.bit_length()-1
            if h in piv: v^=piv[h]
            else: piv[h]=v; r+=1; break
    return r
def span(basis):
    el={0}
    for b in basis: el|={e^b for e in el}
    return el
def enum_sigma(a,b):
    t=a+b; mA=(1<<a)-1; mB=((1<<t)-1)^mA; F=(1<<t)-1
    q=lambda x:(pc(x&mA)-pc(x&mB))%4
    Q=lambda v:(q(v&F)-q(v>>t))%4                     # q (+) -q on V(+)V
    iso=[v for v in range(1,1<<(2*t)) if Q(v)==0]
    one=(1<<(2*t))-1; found=set()
    def dfs(basis,el,start):
        if len(basis)==t: found.add(frozenset(el)); return
        for i in range(start,len(iso)):
            v=iso[i]
            if v in el: continue
            if any(pc(v&u)%2 for u in basis): continue
            dfs(basis+[v], el|{e^v for e in el}, i+1)
    dfs([one],{0,one},0)
    return t,q,[sorted(T) for T in found]
def check(a,b):
    t,q,Sg=enum_sigma(a,b); F=(1<<t)-1; N=len(Sg)
    c=lambda k: prod(2**i+1 for i in range(t-k-1))
    w=min(a,b)+kappa(abs(a-b))
    NL=lambda T:[v&F for v in T if (v>>t)==0]
    NR=lambda T:[v>>t for v in T if (v&F)==0]
    # composition (own): T o T' = {(x,z): exists y, (x,y) in T, (y,z) in T'}
    def comp(T,Tp):
        byy={}
        for u in Tp: byy.setdefault(u&F,[]).append(u>>t)
        R=set()
        for u in T:
            for z in byy.get(u>>t,[]): R.add((u&F)|(z<<t))
        return frozenset(R)
    idx={frozenset(T):i for i,T in enumerate(Sg)}
    D=sp.symbols('D')
    ok_count=(N==c(0)); ok_cols=True; maxlev=0
    for Tp in Sg:
        r=rank(NL(Tp)); maxlev=max(maxlev,r)
        Kp=NL(Tp)
        S=sum(D**(rank(NR(T))+rank(Kp)-rank(NR(T)+Kp)) for T in Sg)   # D^{dim(N_R(T) cap K')}
        ok_cols &= sp.expand(S - c(r)*sp.prod([D+2**i for i in range(t-r-1,t-1)]))==0
    # Lemma 2 on every isotropic K that occurs as some N_L(T')
    ok_L2=True
    for Tp in Sg:
        K=NL(Tp); k=rank(K); el=span(K)
        cnt=sum(1 for T in Sg if el<=set(NR(T)) | {0} and all(x in span(NR(T)) for x in el))
        ok_L2 &= cnt==c(k)
    # spectral radius of M at D=2^n, n=1..3, own composition
    rho_ok=True
    for n in (1,2,3):
        Dv=2**n; Z=Dv*prod(Dv+2**i for i in range(t-1)); M=np.zeros((N,N))
        for j,Tp in enumerate(Sg):
            for T in Sg:
                R=comp(T,Tp); d=rank(NR(T))+rank(NL(Tp))-rank(NR(T)+NL(Tp))
                M[idx[R],j]+=Dv**d/Z
        rho=max(abs(np.linalg.eigvals(M)))
        Fv=prod(2**i+1 for i in range(t-w-1))/(Dv*prod(Dv+2**i for i in range(t-w-1)))
        rho_ok &= abs(rho-Fv)<1e-12
    print(f"(a,b)=({a},{b}) t={t}: |Sigma|={N} (c_0={c(0)}) ok={ok_count}; max level={maxlev} (w={w}); "
          f"column sums = F_(t-r) all columns: {ok_cols}; Lemma 2 counts: {ok_L2}; rho(M)=F_(t-w) at n=1,2,3: {rho_ok}")
for ab in [(2,1),(3,0),(2,2),(3,1),(4,0),(3,2),(4,1),(5,0)]:
    check(*ab)
# Lemma 3 identity, own code, t<=30
D=sp.symbols('D'); bad=0; n=0
qb=lambda r,k: prod(2**(r-i)-1 for i in range(k))//prod(2**(i+1)-1 for i in range(k))
for t in range(2,31):
    c=lambda k: prod(2**i+1 for i in range(t-k-1))
    for r in range(0,t//2+1):
        n+=1; bad+= sp.expand(sum(qb(r,k)*sp.prod([D-2**i for i in range(k)])*c(k) for k in range(r+1)) - c(r)*sp.prod([D+2**i for i in range(t-r-1,t-1)]))!=0
print(f"Lemma 3 q-identity, own code: {n} (t,r) pairs, t<=30, failures {bad}")
