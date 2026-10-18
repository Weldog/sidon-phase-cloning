import numpy as np, pytest
from qhlab import lemma_r as LR, attacks as AT

def cbq(K, q):
    x = np.arange(q); return np.exp(2j * np.pi * np.outer(x, K) / q) / np.sqrt(len(K))

@pytest.mark.parametrize("K,q,T,k", [([0, 1, 5, 17], 64, 1, 1), ([0, 1, 2, 3], 64, 1, 1),
                                     (None, None, 2, 1), (None, None, 1, 2)])
def test_lemma_r_equals_full_sdp(K, q, T, k):
    if K is None:
        K, q = LR.torus_K(2 if (T, k) == (2, 1) else 3, T, k)
    assert abs(LR.lemma_r(K, q, T, k) - AT.counterfeit_sdp(cbq(np.array(K), q), T, k)) < 1e-5

@pytest.mark.parametrize("k", [1, 2, 3])
def test_torus_t1_closed_form_d4(k):
    K, q = LR.torus_K(4, 1, k)
    assert abs(LR.lemma_r(K, q, 1, k) - LR.t1_closed_form(4, k)) < 1e-7

def test_sidon_only_fails_at_t_plus_k_3():
    K, q = [18, 21, 27, 28], 32
    assert LR.is_Bh(K, q, 2) and not LR.is_Bh(K, q, 3)
    assert abs(LR.lemma_r(K, q, 1, 1) - 7 / 16) < 1e-7          # Sidon => torus at (1,1)
    assert LR.lemma_r(K, q, 1, 2) > LR.t1_closed_form(4, 2) + 0.01  # but leaks at (1,2)

def test_theorem1_kraus_channel_any_K():
    for K, q in [([0, 1, 5, 17], 64), ([0, 1, 2, 3], 64)]:
        assert abs(LR.sidon_kraus_fidelity(cbq(np.array(K), q)) - 7 / 16) < 1e-12

def test_bose_chowla_p7_h3():
    K, q = LR.bose_chowla(7, 3)
    assert q == 342 and LR.is_Bh(K, q, 3) and not LR.is_Bh(K, q, 4)
    ind = np.zeros(q); np.add.at(ind, K, 1)
    assert abs((np.abs(np.fft.fft(ind))[1:].max() / 7) ** 2 - 0.525011) < 1e-5
    assert abs(LR.lemma_r(K, q, 1, 2) - LR.t1_closed_form(7, 2)) < 1e-7   # B_3 => torus at (1,2)
