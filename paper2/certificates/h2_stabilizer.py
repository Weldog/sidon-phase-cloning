# Independent shim (Claude): occupation lists and stabilizer states, written from scratch (not Grok's module).
import itertools, numpy as np
def occs(d, n):
    return [tuple(c.count(i) for i in range(d)) for c in itertools.combinations_with_replacement(range(d), n)]
def stab_states(nq):
    H=np.array([[1,1],[1,-1]])/np.sqrt(2); S=np.diag([1,1j]); I2=np.eye(2)
    def on(q,g):
        m=np.array([[1]],complex)
        for j in range(nq): m=np.kron(m, g if j==q else I2)
        return m
    gens=[on(q,H) for q in range(nq)]+[on(q,S) for q in range(nq)]
    D=2**nq
    for a in range(nq):
        for b in range(nq):
            if a!=b:
                C=np.zeros((D,D))
                for x in range(D):
                    bits=[(x>>(nq-1-j))&1 for j in range(nq)]
                    if bits[a]: bits[b]^=1
                    y=sum(bit<<(nq-1-j) for j,bit in enumerate(bits)); C[y,x]=1
                gens.append(C)
    def canon(v):
        i=np.argmax(np.abs(v)>1e-9); v=v*np.exp(-1j*np.angle(v[i])); return tuple(np.round(v,8))
    st=np.zeros(D,complex); st[0]=1; seen={canon(st):st}; fr=[st]
    while fr:
        nf=[]
        for v in fr:
            for g in gens:
                w=g@v; c=canon(w)
                if c not in seen: seen[c]=w; nf.append(w)
        fr=nf
    return list(seen.values())
