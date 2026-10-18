import numpy as np, importlib.util, pathlib, itertools
src = (pathlib.Path(__file__).parent.parent / "r6_verify.py").read_text().split("out = {}")[0]
ns = {}; exec(src, ns)

def test_theorem_F_fusion_equals_additive_energy():
    exec((pathlib.Path(__file__).parent.parent / "r6_verify.py").read_text().split("# ---------- (1) Theorem F ----------")[1]
         .split('out["thmF"]')[0], ns)
    v, f = ns["fusion"]([0, 1, 3], 7)
    assert abs(v - f) < 1e-6 and abs(f - 5 / 9) < 1e-12

def _real_states(P): return np.asarray(P, complex)

def test_conjecture_A_noncoherent_icosa_11_and_hierarchy_counterexample():
    phi = (1 + 5 ** 0.5) / 2
    ax = [np.array(v, float) / np.linalg.norm(v) for v in
          [(0, 1, phi), (0, 1, -phi), (1, phi, 0), (1, -phi, 0), (phi, 0, 1), (-phi, 0, 1)]]   # 6 icosahedral axes
    w = np.ones(6) / 6
    z, wz = np.polynomial.legendre.leggauss(6); Q, W = [], []
    for zi, wi in zip(z, wz):
        for j in range(12):
            p = 2 * np.pi * j / 12; r = np.sqrt(1 - zi * zi); Q.append([r * np.cos(p), r * np.sin(p), zi]); W.append(wi / 24)
    PH11 = ns["counterfeit"](_real_states(Q), np.array(W), 1, 1)
    PX11 = ns["counterfeit"](_real_states(ax), w, 1, 1)
    assert abs(PX11 - PH11) < 1e-7 and abs(PH11 - 0.6) < 1e-7               # moment-matched, mixed not: equality (Conj. A 'if')
    PX12 = ns["counterfeit"](_real_states(ax), w, 1, 2)
    assert abs(PX12 - 0.6) < 1e-7                                            # NOT moment-matched yet = P_H: hierarchy 'only if' false
