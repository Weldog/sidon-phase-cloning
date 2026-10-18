"""Independent audit of Grok Round 4 using QH-Lab code (not Grok's)."""
import json, numpy as np
from qhlab import designs as Ds, stabilizer as ST, attacks as AT, theorems as TH, phase_hash as PH
AT.SOLVER.update(eps=1e-7)
out = {}
s1, s2, s3 = (ST.stabilizer_states(n) for n in (1, 2, 3))
# 1. S-W break, n=1 by full Choi SDP (our code) vs Grok table vs Theorem A
rows = []
for (T, k, grok) in [(1, 3, 0.416667), (2, 2, 0.625), (3, 1, 0.833333)]:
    rows.append(dict(n=1, T=T, k=k, grok=grok, sdp=AT.counterfeit_sdp(s1, T, k), thmA=TH.theorem_a(s1, T, k),
                     W=Ds.werner(2, T, k)))
# 2. Theorem A at n=2,3 vs Grok table
for (S, n, d, T, k, grok) in [(s2, 2, 4, 1, 3, 0.125), (s2, 2, 4, 2, 2, 0.3125), (s2, 2, 4, 3, 1, 0.625),
                              (s3, 3, 8, 1, 3, 0.027778), (s3, 3, 8, 3, 1, 0.416667)]:
    rows.append(dict(n=n, T=T, k=k, grok=grok, thmA=TH.theorem_a(S, T, k), W=Ds.werner(d, T, k)))
out["SW_break"] = rows
# 3. Theorem B certificate: exact W for T+k<=3 on stabilizer states
out["thmB_stab_n2"] = {f"{T},{k}": (TH.certificate_b(s2, T, k), Ds.werner(4, T, k)) for T, k in [(1, 1), (1, 2), (2, 1), (1, 3)]}
# 4. Witting-40
V = TH.witting40()
out["witting"] = None if V is None else dict(N=len(V), delta2=PH.delta_bruteforce(V) ** 2,
    t3=Ds.is_t_design(V, 3, 1e-9), t4=Ds.is_t_design(V, 4, 1e-9),
    cf_1to2_sdp=AT.counterfeit_sdp(V, 1, 1), cf_1to3_thmB=TH.certificate_b(V, 1, 2), W12=Ds.werner(4, 1, 1), W13=Ds.werner(4, 1, 2),
    levenstein=PH.levenstein_floor(40, 4))
# 5. Phase codes: Grok's generic formula (2d-1)/d^2 and the Sidon hypothesis
ph = []
for (L, K) in [(6, [0, 1, 5, 17]), (6, [10, 30, 45, 60]), (6, [0, 1, 2, 3]), (7, [0, 1, 4, 9, 11]), (7, [0, 1, 4, 10, 12, 17])]:
    S = PH.codebook(L, K); d = len(K)
    ph.append(dict(L=L, K=K, d=d, sidon=TH.is_sidon(K, 2 ** L), sdp_1to2=AT.counterfeit_sdp(S, 1, 1),
                   formula=(2 * d - 1) / d ** 2, W=Ds.werner(d, 1, 1), delta2=PH.delta_fft(L, K) ** 2))
out["phase"] = ph
json.dump(out, open("r4_verify.json", "w"), indent=1, default=float)
for k, v in out.items():
    print(k); print(json.dumps(v, indent=0, default=float)[:2500])
