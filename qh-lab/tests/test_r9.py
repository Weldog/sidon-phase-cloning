import numpy as np, importlib.util, pathlib
from qhlab import lemma_r as LR
spec = importlib.util.spec_from_file_location("r9", pathlib.Path(__file__).parent.parent / "r9_verify.py")
r9 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r9)

def test_spin1_octahedral_kill_grok_fiducial():
    f = np.array([0.535964010446345 - 0.638498468033802j, 0.402234308639641 - 0.222141471386666j,
                  -0.298671263020807 - 0.068691266647742j]); f /= np.linalg.norm(f)
    r, RX, RH = r9.evaluate(f)
    assert r["M2_diff"] < 1e-12 and r["thmB"] > r["P_H"] + 0.1
    assert r["P_X"] - r["P_H"] > 6.4e-4                     # Grok: certified >= 6.437e-4

def test_torus_minimal_q_B_h_codebook_equals_torus_d2_large_T():
    K, q = [0, 1], 14                                        # B_13 but not B_25 mod 14
    assert LR.is_Bh(K, q, 13) and not LR.is_Bh(K, q, 25)
    Kt, qt = LR.torus_K(2, 12, 1)
    assert abs(LR.lemma_r(K, q, 12, 1) - LR.lemma_r(Kt, qt, 12, 1)) < 1e-7
