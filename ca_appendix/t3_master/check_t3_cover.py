#!/usr/bin/env python3
"""R11: closes Theorem T3. Checks that (1) t3_bounds.py reported all analytic regions < 1 and (2) every finite pair
(d,k) with S3(d,k) >= 1 has an exact clique-LP certificate at T = 3 in ../t3_clique_lp/cert.json whose verified
bound (re-run ../t3_clique_lp/check_torus_lp.py; expected_output.txt there) is <= 1. Pure bookkeeping, no arithmetic."""
import json, re, sys
b = json.load(open('t3_bounds_DU100.json'))
cert = {(c['d'], c['T'], c['k']) for c in json.load(open('../t3_clique_lp/cert.json'))}
passed = {}
for line in open('../t3_clique_lp/expected_output.txt'):
    m = re.match(r'd=(\d+) T=(\d+) k=(\d+)\s+bound<=([0-9.]+).*PASS', line)
    if m: passed[(int(m[1]), int(m[2]), int(m[3]))] = float(m[4])
ok = b['ANALYTIC_REGIONS_OK'] and b['E3_max_at_m0'] < 1 and b['S3_max_below_1'][0] < 1
res = [tuple(x) for x in b['S3_ge_1']]
for d, k in res:
    key = (d, 3, k); good = key in cert and key in passed and passed[key] <= 1; ok &= good
    print('residual (d,k)=(%d,%d): certificate %s, verified bound %s' % (d, k, key in cert, passed.get(key)))
print('residual pairs: %d ; max verified LP bound on them: %.10f' % (len(res), max(passed[(d, 3, k)] for d, k in res)))
print('THEOREM T3 (all d>=3, k>=1) COMPUTATIONAL STEPS:', 'PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
