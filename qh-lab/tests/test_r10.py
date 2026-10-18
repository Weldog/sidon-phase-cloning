import numpy as np, importlib.util, pathlib
from scipy.optimize import least_squares
from qhlab import lemma_r as LR

def test_T2_torus_theorem_d12_k1_B3_not_B5():
    d, T, k = 12, 2, 1; K = [4 ** j for j in range(d)]; q = 3 * 4 ** (d - 1) + 1
    Kt, qt = LR.torus_K(d, T, k)
    assert abs(LR.lemma_r(K, q, T, k) - LR.lemma_r(Kt, qt, T, k)) < 1e-8

def test_spin52_binary_icosahedral_window_prediction_fail():
    spec = importlib.util.spec_from_file_location("s52", pathlib.Path(__file__).parent.parent / "r10_spin52.py")
    s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
    assert abs(s.irr - 1) < 1e-12
    x0 = np.random.default_rng(900).normal(size=12)
    f_of = lambda x: (lambda p: p / np.linalg.norm(p))(x[:6] + 1j * x[6:])
    def resid(x):
        f = f_of(x); Dm = s.M2(np.array([U @ f for U in s.DG]), s.w60) - s.M2(s.Uq @ f, s.Wq)
        return np.concatenate([Dm.real.ravel(), Dm.imag.ravel()])
    f = f_of(least_squares(resid, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000).x)
    X = np.array([U @ f for U in s.DG])
    assert np.abs(s.M2(X, s.w60) - s.M2(s.Uq @ f, s.Wq)).max() < 1e-12
    assert 6 * s.lmax(s.Rm_(X, s.w60)) - 6 * s.lmax(s.Rm_(s.Uq @ f, s.Wq)) > 1e-4
