import numpy as np, importlib.util, pathlib
spec = importlib.util.spec_from_file_location("r7", pathlib.Path(__file__).parent.parent / "r7_verify.py")
r7 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r7)

def test_spin32_ghz_genuine_pass_and_untuned_control():
    fid = np.array([1, 0, 0, 1], complex) / np.sqrt(2)
    H = np.array([U @ fid for U in r7.Uq]); MH = r7.M(H, r7.Wq, 2); w60 = np.ones(60) / 60
    orbit = lambda e3: np.array([Dg @ r7.D_euler(*e3) @ fid for Dg in r7.DG])
    def resid(e3):
        Dm = r7.M(orbit(e3), w60, 2) - MH; return np.concatenate([Dm.real.ravel(), Dm.imag.ravel()])
    x = r7.least_squares(resid, (0.1, 0.4, 0.2), xtol=1e-15, ftol=1e-15, gtol=1e-15).x
    X = orbit(x)
    assert np.abs(r7.M(X, w60, 2) - MH).max() < 1e-10                       # (T+k)-moment matched
    assert np.abs(r7.Rmix(X, w60, 1, 1) - r7.Rmix(H, r7.Wq, 1, 1)).max() > 1e-3  # mixed moment NOT matched
    PH = r7.P_opt(H, r7.Wq, 1, 1)
    assert abs(PH - 31 / 70) < 1e-7 and r7.refinedB(X, w60, 1, 1) > PH + 0.1   # Theorem B loose -> genuine
    assert abs(r7.P_opt(X, w60, 1, 1) - PH) < 1e-7                            # Conjecture A holds
    assert r7.P_opt(orbit((0.3, 0.7, 1.1)), w60, 1, 1) > PH + 1e-3           # untuned control fails strictly
