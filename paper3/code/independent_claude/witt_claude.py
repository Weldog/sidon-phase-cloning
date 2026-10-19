# Claude's independent check of the Witt index of q_{a,b}(x) = |x_A| - |x_B| mod 4 on F_2^t (own code, not Grok's).
import sys, time
def pc(v): return bin(v).count('1')
def maxiso(a,b,cap):
    t=a+b; mA=(1<<a)-1; mB=((1<<t)-1)^mA
    q=lambda x:(pc(x&mA)-pc(x&mB))%4
    iso=[v for v in range(1,1<<t) if q(v)==0]
    best=[0]
    def span_add(basis,v):
        # reduce v against basis (dict pivot->vec); return reduced
        for p,bv in basis.items():
            if v>>p &1: v^=bv
        return v
    def dfs(basis,cands,depth):
        if depth>best[0]: best[0]=depth
        if best[0]>=cap: return
        if depth+len(cands)<=best[0]: return
        for i,v in enumerate(cands):
            r=span_add(basis,v)
            if r==0: continue
            p=r.bit_length()-1
            nb=dict(basis); 
            for pp in list(nb):
                if nb[pp]>>p &1: nb[pp]^=r
            nb[p]=r
            nc=[c for c in cands[i+1:] if pc(c&v)%2==0]
            dfs(nb,nc,depth+1)
            if best[0]>=cap: return
    dfs({},iso,0); return best[0]
mu=[0,1,2,3,2,3,2,1]
bad=0; n=0; t0=time.time()
for t in range(1,10):
    for b in range(0,t+1):
        a=t-b; pred=(t-mu[(a-b)%8])//2
        got=maxiso(a,b,pred+1)   # search up to pred+1: finding pred+1 would refute the upper bound
        n+=1; bad+= got!=pred
        if got!=pred: print("MISMATCH",a,b,got,pred)
print(f"Witt index, all (a,b) with t=a+b<=9: {n} forms, mismatches {bad}  [{time.time()-t0:.0f}s]")
