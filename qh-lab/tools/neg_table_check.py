"""Reproduce every value of the paper's negative-results table with the QH-Lab reduced dual (Lemma R).
Usage: python3 tools/neg_table_check.py   (runtime ~1-4 min; prints max deviation and the per-row worst set)."""
import json, os, sys, collections
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from qhlab import lemma_r as LR
D = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'data', 'neg_table.json')))
def find(d, T, k, typ):
    for s in D['sets']:
        if s['d'] == d and s['type'] == typ and s['tk'] in (f"({T},{k})", 'all'): return s['K'], s['q']
mx, worst = 0.0, collections.Counter()
for r in D['rows']:
    d, T, k = r['d'], r['T'], r['k']; vals = {}
    Kt, qt = LR.torus_K(d, T, k); vals['tor'] = LR.lemma_r(Kt, qt, T, k)
    for typ in ('sidon', 'singer', 'ap', 'random'):
        K, q = find(d, T, k, typ); vals[typ] = LR.lemma_r(K, q, T, k)
    for key, v in vals.items(): mx = max(mx, abs(v - r[key]))
    worst[max(('sidon', 'singer', 'ap', 'random'), key=lambda x: vals[x])] += 1
    assert vals['sidon'] > vals['tor'] and vals['singer'] > vals['tor']
print(f"{len(D['rows'])} rows x 5 values; max |QH-Lab - table| = {mx:.2e}; worst set per row: {dict(worst)}")
assert mx < 1e-6
