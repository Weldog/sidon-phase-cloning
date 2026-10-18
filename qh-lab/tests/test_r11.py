import pytest
from fractions import Fraction
from qhlab import lemma_r as LR

def test_T3_torus_theorem_spot_d4_k1():
    d, T, k = 4, 3, 1; K = [5 ** j for j in range(d)]; q = 4 * 5 ** (d - 1) + 1
    Kt, qt = LR.torus_K(d, T, k)
    assert abs(LR.lemma_r(K, q, T, k) - LR.lemma_r(Kt, qt, T, k)) < 1e-8

def test_T1_k_gt_d_matches_eq16_and_closed_form_guarded():
    Kt, qt = LR.torus_K(2, 1, 3)
    assert abs(LR.lemma_r(Kt, qt, 1, 3) - 5 / 8) < 1e-8          # arXiv:2610.06723 Eq. (16): 5/8
    with pytest.raises(ValueError): LR.t1_closed_form(2, 3)
