"""T=2 Torus Theorem spot checks (independent QH-Lab sparse Lemma R): B_{2+k} codebooks that are NOT B_{4+k}, vs torus."""
import time, json
from qhlab import lemma_r as LR
out = []
for d, k in [(5, 1), (5, 3), (8, 1), (8, 2), (12, 1), (12, 2), (16, 1), (20, 1)]:
    T = 2; h = T + k; b = h + 1
    K = [b ** j for j in range(d)]; q = h * b ** (d - 1) + 1
    if d <= 8: assert LR.is_Bh(K, q, h) and not LR.is_Bh(K, q, 2 * T + k)
    t = time.time(); pc = LR.lemma_r(K, q, T, k); Kt, qt = LR.torus_K(d, T, k); pt = LR.lemma_r(Kt, qt, T, k)
    out.append(dict(d=d, T=T, k=k, P_code=pc, P_torus=pt, diff=pc - pt, sec=round(time.time() - t, 1))); print(out[-1], flush=True)
json.dump(out, open("t2_spot.json", "w"), indent=1)
