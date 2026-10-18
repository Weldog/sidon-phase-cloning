"""B1 run: reproduce Round-3 numbers, then evaluate QDH-2 arms (stabilizer, Clifford-orbit) vs QDH-1 phase hash."""
import json, time
import numpy as np
from qhlab import designs as Ds, phase_hash as PH, stabilizer as ST, attacks as AT, clifford as CL

rng = np.random.default_rng(20261015)
R = {"meta": {"version": "0.1.0-b1", "seed": 20261015, "date": "2026-10-15"}}
t0 = time.time()

def eq(N): return np.array([[1, np.exp(2j*np.pi*x/N)] for x in range(N)])/np.sqrt(2)
a = np.arccos(-1/3)
SIC2 = np.array([[1,0]]+[[np.cos(a/2), np.exp(2j*np.pi*k/3)*np.sin(a/2)] for k in range(3)], complex)

# A. Reproduction of Round-3 (Grok) numbers
R["A_reproduction"] = [
  {"claim":"SIC d2 preimage T=2","grok":0.750,"ours":AT.preimage_sdp(SIC2,2)},
  {"claim":"SIC d2 preimage T=3","grok":0.9714,"ours":AT.preimage_sdp(SIC2,3)},
  {"claim":"Phase d2 N8 preimage T=2","grok":0.3643,"ours":AT.preimage_sdp(eq(8),2)},
  {"claim":"Phase d2 N8 counterfeit 1->2","grok":0.750,"ours":AT.counterfeit_sdp(eq(8),1,1)},
  {"claim":"Phase d2 N16 counterfeit 1->2","grok":0.750,"ours":AT.counterfeit_sdp(eq(16),1,1)},
  {"claim":"Phase d2 N16 counterfeit 1->3","grok":0.750,"ours":AT.counterfeit_sdp(eq(16),1,2)},
  {"claim":"Phase d2 N8 eta_{1,1}","grok":0.1405,"ours":Ds.design_defect(eq(8),1,1)},
]

# B. Stabilizer arm (exact 3-design)
stab = {}
for n in (1, 2, 3):
    S = ST.stabilizer_states(n); d = 2**n; N = len(S)
    row = {"n":n,"d":d,"N":N,"L_bits":float(np.log2(N)),
           "t3_design":Ds.is_t_design(S,3),"t4_design":Ds.is_t_design(S,4),
           "delta2":PH.delta_bruteforce(S)**2,
           "preimage_pgm":{T:AT.pgm_success(S,T) for T in (1,2,3,4)},
           "preimage_bound":{T:min(1,Ds.D(d,T)/N) for T in (1,2,3,4)},
           "werner_1to2":Ds.werner(d,1,1)}
    if n <= 2:
        row["counterfeit_1to2_sdp"] = AT.counterfeit_sdp(S,1,1)
    if n == 2:
        row["counterfeit_1to3_sdp"] = AT.counterfeit_sdp(S,1,2); row["werner_1to3"] = Ds.werner(4,1,2)
        row["eta_1_2"] = Ds.design_defect(S,1,2)
    else:
        row["counterfeit_1to2_note"] = "eta_{1,1}=0 (exact 3-design) => Thm part 2 gives P = W exactly; SDP not needed"
    if n == 1:
        row["counterfeit_1to3_sdp"] = AT.counterfeit_sdp(S,1,2); row["werner_1to3"] = Ds.werner(2,1,2)
        row["counterfeit_2to3_sdp"] = AT.counterfeit_sdp(S,2,1); row["werner_2to3"] = Ds.werner(2,2,1)
        row["eta_1_2"] = Ds.design_defect(S,1,2); row["eta_2_1"] = Ds.design_defect(S,2,1)
    stab[n] = row
R["B_stabilizer"] = stab

# C. Phase hash at matched dimension d=4
ph = []
for L in (6, 13):
    best = None
    for _ in range(200):
        K = PH.random_K(L, 4, rng); dl = PH.delta_fft(L, K)
        if best is None or dl < best[0]: best = (dl, K)
    S = PH.codebook(L, best[1])
    ph.append({"L":L,"N":2**L,"d":4,"K":best[1].tolist(),"delta2":best[0]**2,
               "eta_1_1":Ds.design_defect(S,1,1),"counterfeit_1to2_sdp":AT.counterfeit_sdp(S,1,1),
               "werner_1to2":Ds.werner(4,1,1),"welch":PH.welch_floor(2**L,4),
               "levenstein":PH.levenstein_floor(2**L,4)})
R["C_phase_d4"] = ph

# D. Clifford-orbit arm: 3-design for ANY fiducial; tune fiducial to minimise collision bias delta
co = {}
for n in (1, 2):
    G = CL.clifford_group(n); d = 2**n
    dl, f = CL.best_fiducial(G, trials=3000 if n == 1 else 1500, rng=rng)
    S = CL.orbit_states(G, f); N = len(S)
    row = {"n":n,"d":d,"N":N,"delta2":dl**2,"t3_design":Ds.is_t_design(S,3,1e-8),
           "welch":PH.welch_floor(N,d),"levenstein":PH.levenstein_floor(N,d),
           "werner_1to2":Ds.werner(d,1,1),"fiducial":[[c.real,c.imag] for c in f]}
    if n == 1:
        row["counterfeit_1to2_sdp"] = AT.counterfeit_sdp(S,1,1)
        row["eta_1_1"] = Ds.design_defect(S,1,1)
    co[n] = row
R["D_clifford_orbit"] = co
R["meta"]["runtime_s"] = round(time.time()-t0,1)
json.dump(R, open("results.json","w"), indent=1, default=lambda o: o.item() if hasattr(o,"item") else str(o))
print(json.dumps(R, indent=1, default=lambda o: o.item() if hasattr(o,"item") else str(o))[:6000])
