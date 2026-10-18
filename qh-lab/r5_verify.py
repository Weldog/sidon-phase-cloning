import json, numpy as np
from qhlab import lemma_r as LR, attacks as AT, designs as Ds, phase_hash as PH, theorems as TH
AT.SOLVER.update(eps=1e-8)
out = {}

def codebook_q(K, q):
    x = np.arange(q); return np.exp(2j*np.pi*np.outer(x, K)/q)/np.sqrt(len(K))
# 1. Lemma R (ours) vs full Choi SDP (ours), across (T,k)
cases = [("Sidon d4", [0,1,5,17], 64, 1, 1), ("AP d4", [0,1,2,3], 64, 1, 1), ("coincid d4", [10,30,45,60], 64, 1, 1)]
for d, T, k in [(2,2,1), (2,1,3), (2,3,1), (3,1,2), (4,1,2)]:
    K, q = LR.torus_K(d, T, k); cases.append((f"torus d{d}", K.tolist(), q, T, k))
rows = []
for name, K, q, T, k in cases:
    S = codebook_q(np.array(K), q)
    rows.append(dict(case=name, K=K, q=q, T=T, k=k, lemmaR=LR.lemma_r(K, q, T, k), full_sdp=AT.counterfeit_sdp(S, T, k)))
out["lemmaR_vs_fullSDP"] = rows
# 2. d=4 torus values via Lemma R vs Grok table
grok = {(1,1):0.4375,(1,2):0.28125,(2,1):0.577816,(1,3):0.234375,(2,2):0.428776,(3,1):0.682631}
out["torus_d4"] = [dict(T=T, k=k, ours=LR.lemma_r(*LR.torus_K(4,T,k), T, k), grok=g,
                        T1_formula=(LR.t1_closed_form(4,k) if T==1 else None), W=Ds.werner(4,T,k),
                        thm2_bound=Ds.werner(4,T,k)*np.prod([1+i/4 for i in range(T+k)])) for (T,k), g in grok.items()]
# 3. Sidon-only counterexamples (B2, not B3) at d=4, q=32
rng = np.random.default_rng(7); ce = []
while len(ce) < 3:
    K = sorted(rng.choice(32, 4, replace=False).tolist())
    if LR.is_Bh(K, 32, 2) and not LR.is_Bh(K, 32, 3):
        r = dict(K=K, q=32, P12=LR.lemma_r(K,32,1,2), P21=LR.lemma_r(K,32,2,1), P11=LR.lemma_r(K,32,1,1))
        if len(ce) == 0: r["P12_fullSDP"] = AT.counterfeit_sdp(codebook_q(np.array(K),32),1,2)
        ce.append(r)
out["sidon_only_d4"] = dict(torus12=0.28125, torus21=0.577816, torus11=0.4375, sets=ce)
# 4. Theorem-1 primal channel attains (2d-1)/d^2 on every phase codebook (any K)
out["thm1_kraus"] = [dict(K=K, q=q, fid=LR.sidon_kraus_fidelity(codebook_q(np.array(K), q)), target=(2*len(K)-1)/len(K)**2)
                     for K, q in [([0,1,5,17],64), ([0,1,2,3],64), ([0,1,4,9,11],128)]]
# 5. Theorem 2 certificate on a B_3 set (= W*(1+1/d)(1+2/d)) and Bose-Chowla
out["bose_chowla"] = []
for p, h in [(7,3), (13,3)]:
    K, q = LR.bose_chowla(p, h)
    out["bose_chowla"].append(dict(p=p, h=h, q=q, d=len(K), B2=LR.is_Bh(K,q,2), B3=LR.is_Bh(K,q,3), B4=LR.is_Bh(K,q,4),
                                   delta2=PH.delta_fft(int(np.log2(q)) if False else 0, []) if False else None))
for r in out["bose_chowla"]:
    K, q = LR.bose_chowla(r["p"], r["h"])
    ind = np.zeros(q); np.add.at(ind, K, 1); f = np.abs(np.fft.fft(ind))[1:] / len(K)
    r["delta2"] = float(f.max() ** 2); r["katz"] = (r["h"] - 1) ** 2 / r["p"]
    if r["p"] == 7:
        r["thm2_cert_12"] = TH.certificate_b(codebook_q(K, q), 1, 2)
        r["thm2_formula_12"] = Ds.werner(7,1,2) * (1+1/7) * (1+2/7)
        r["P12_lemmaR"] = LR.lemma_r(K, q, 1, 2); r["torus12_T1formula"] = LR.t1_closed_form(7, 2)
json.dump(out, open("r5_verify.json","w"), indent=1, default=float)
