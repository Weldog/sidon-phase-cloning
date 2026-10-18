import numpy as np, importlib.util, pathlib
from scipy.optimize import least_squares
spec = importlib.util.spec_from_file_location("r8", pathlib.Path(__file__).parent.parent / "r8_verify.py")
r8 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r8)

def test_su3_sym2_hessian_counterexample_to_conjecture_A():
    assert len(r8.HG) == 216 and abs(r8.irr - 1) < 1e-12                 # Hessian group, irreducible on Sym^2 C^3
    x0 = np.random.default_rng(100).normal(size=12)
    phi_of = lambda x: (lambda p: p / np.linalg.norm(p))(x[:6] + 1j * x[6:])
    def resid(x):
        Dm = r8.M2_X(phi_of(x)) - r8.M2_H(phi_of(x)); return np.concatenate([Dm.real.ravel(), Dm.imag.ravel()])
    phi = phi_of(least_squares(resid, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000).x)
    M2X = r8.M2_X(phi)
    assert np.abs(M2X - r8.M2_H(phi)).max() < 1e-12                       # hypothesis of Conjecture A holds
    PX, PH = 6 * r8.lmax(r8.R_X(phi)), 6 * r8.lmax(r8.R_H(phi))
    assert 6 * r8.lmax(M2X) > PH + 0.05                                    # Theorem B loose -> genuine
    assert PX - PH > 1e-6                                                  # ...and the conclusion fails
