import numpy as np, pytest
from qhlab import designs, phase_hash, stabilizer, attacks, fib, delaunay

def eq(N):
    return np.array([[1, np.exp(2j * np.pi * x / N)] for x in range(N)]) / np.sqrt(2)

def sic2():
    a = np.arccos(-1 / 3)
    S = [[1, 0]] + [[np.cos(a / 2), np.exp(2j * np.pi * k / 3) * np.sin(a / 2)] for k in range(3)]
    return np.array(S, complex)

# --- T1: Fibonacci representation consistency (benchmark module) ---
@pytest.mark.parametrize("n,total,dim", [(3, 1, 2), (4, 0, 2), (6, 0, 5), (7, 1, 13), (9, 1, 34)])
def test_fib_dims_unitary_braid_relations(n, total, dim):
    g = fib.generators(n, total)
    assert g[0].shape[0] == dim
    I = np.eye(dim)
    for s in g:
        assert np.allclose(s @ s.conj().T, I, atol=1e-12)
        assert np.allclose(np.linalg.matrix_power(s, 10), I, atol=1e-10)   # R^10 = 1
    for i in range(len(g) - 1):
        a, b = g[i], g[i + 1]
        assert np.allclose(a @ b @ a, b @ a @ b, atol=1e-12)               # Yang-Baxter
    for i in range(len(g)):
        for j in range(i + 2, len(g)):
            assert np.allclose(g[i] @ g[j], g[j] @ g[i], atol=1e-12)        # far commutation

# --- phase hash: FFT certificate == brute force ---
@pytest.mark.parametrize("L,d,seed", [(6, 3, 1), (8, 5, 2), (10, 8, 3)])
def test_fft_delta_matches_bruteforce(L, d, seed):
    K = phase_hash.random_K(L, d, np.random.default_rng(seed))
    assert abs(phase_hash.delta_fft(L, K) - phase_hash.delta_bruteforce(phase_hash.codebook(L, K))) < 1e-9

# --- designs ---
def test_stabilizer_counts_and_design_order():
    s1, s2 = stabilizer.stabilizer_states(1), stabilizer.stabilizer_states(2)
    assert len(s1) == 6 and len(s2) == 60
    for s in (s1, s2):
        assert designs.is_t_design(s, 3) and not designs.is_t_design(s, 4)

def test_sic_is_2_not_3_design():
    assert designs.is_t_design(sic2(), 2) and not designs.is_t_design(sic2(), 3)

def test_haar_moment_trace():
    # Tr X_Haar must equal 1 (average of unit-norm rank-1 projectors)
    assert abs(np.trace(designs.haar_moment_operator(2, 1, 1)) - 1) < 1e-12

def test_design_defect_zero_for_3design_at_m3():
    assert designs.design_defect(stabilizer.stabilizer_states(1), 1, 1) < 1e-12

# --- Round-3 reproduction (Grok numbers) ---
def test_reproduce_grok_preimage():
    assert abs(attacks.preimage_sdp(sic2(), 2) - 0.75) < 1e-4
    assert abs(attacks.preimage_sdp(sic2(), 3) - 0.9714) < 1e-3
    assert abs(attacks.pgm_success(eq(8), 2) - 0.3643) < 1e-3

def test_reproduce_grok_counterfeit_and_eta():
    assert abs(attacks.counterfeit_sdp(eq(8), 1, 1) - 0.75) < 1e-4
    assert abs(designs.design_defect(eq(8), 1, 1) - 0.1405) < 1e-3

def test_stabilizer_counterfeit_hits_werner():
    assert abs(attacks.counterfeit_sdp(stabilizer.stabilizer_states(1), 1, 1) - designs.werner(2, 1, 1)) < 1e-4

# --- Delaunay d=2 unit test ---
def test_delaunay_closest_pair_edge():
    rng = np.random.default_rng(0)
    for _ in range(100):
        v = rng.normal(size=(64, 2)) + 1j * rng.normal(size=(64, 2))
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        assert delaunay.closest_pair_is_delaunay_edge(v)
