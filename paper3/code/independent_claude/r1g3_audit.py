import sys, numpy as np
from fractions import Fraction as Fr
sys.path.insert(0,'/home/claude/r22'); import h2_stabilizer as h2
def kappa(m): q,r=divmod(m,8); return 4*q+[0,0,0,0,1,1,2,3][r]
def F(n,t):
    num=Fr(1); den=Fr(2**n)
    for i in range(t-1): num*=(2**i+1); den*=(2**n+2**i)
    return num/den
def Nn(n):
    r=2**n
    for i in range(1,n+1): r*=(2**i+1)
    return r
# (c) exact limit identity: for t >= n+1, F_t = (1/N) * prod_{i=t-1-n}^{t-2} (1 + 2^-i)
bad=0
for n in range(1,7):
    for t in range(n+1,40):
        rhs=Fr(1,Nn(n))
        for i in range(t-1-n,t-1): rhs*= (1+Fr(1,2**i))
        bad+= F(n,t)!=rhs
print("(c) limit identity, n=1..6, t up to 39: mismatches =",bad)
S={n:np.array(h2.stab_states(n)) for n in (1,2,3)}
O={n:S[n].conj()@S[n].T for n in S}
def spec(n,T,k):
    G=(O[n]**k)*(np.abs(O[n])**(2*T)); return np.sort(np.linalg.eigvalsh((G+G.conj().T)/2))[::-1]/len(S[n])
# (a) T=0 formula
bad=0; cnt=0
for n in (1,2,3):
    for k in range(1,13):
        lam=spec(n,0,k)[0]; cnt+=1; bad+= abs(lam-float(F(n,k-kappa(k))))>1e-10
print(f"(a) T=0: lambda_max/N = F_(k-kappa(k)) on {cnt} cases (n=1..3, k=1..12): mismatches = {bad}")
# (b) claimed unitary equivalence for T+k = 0 mod 4: mixed (T,k) spectrum == unmixed order 2T+k spectrum
for n in (2,3):
    for (T,k) in [(1,3),(2,2),(3,1),(1,7),(2,6),(3,5),(1,4),(2,3)]:
        a=spec(n,T,k); b=spec(n,0,2*T+k); same=np.allclose(a,b,atol=1e-10)
        print(f"(b) n={n} (T,k)=({T},{k}) T+k={T+k} mod4={(T+k)%4}: full spectrum equal to unmixed order {2*T+k}? {same}")
