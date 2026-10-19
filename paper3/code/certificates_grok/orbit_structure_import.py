import sys, random
sys.path.insert(0,'/workspace/r25')
from lagr import *
def orbits(a, b, seed=1):
    S = Space(a, b); t = S.t; Sg = sigma(S); idx = {T: i for i, T in enumerate(Sg)}
    gens = []
    def perm_img(p): return [1 << p[i] for i in range(t)]
    for i in range(a-1): p = list(range(t)); p[i], p[i+1] = p[i+1], p[i]; gens.append(graph(S, perm_img(p)))
    if a >= 2: p = list(range(a))[1:] + [0] + list(range(a, t)); gens.append(graph(S, perm_img(p)))
    for i in range(a, t-1): p = list(range(t)); p[i], p[i+1] = p[i+1], p[i]; gens.append(graph(S, perm_img(p)))
    inv = [T for T in Sg if len(left_defect(S, T)) == 1 and len(left_defect(S, tuple(((u & S.full) << t) | (u >> t) for u in T))) == 1]
    rnd = random.Random(seed); gens += rnd.sample(inv, min(4, len(inv)))
    par = list(range(len(Sg)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for i, T in enumerate(Sg):
        for g in gens:
            for R, d in (compose(S, g, T), compose(S, T, g)):
                ri, rj = find(i), find(idx[R])
                if ri != rj: par[ri] = rj
    roots = sorted(set(find(i) for i in range(len(Sg))))
    return S, [Sg[r] for r in roots]

def orbit_levels(a,b):
    S, reps = orbits(a,b)
    return [len(left_defect(S,T)).bit_length()-1 for T in reps]
