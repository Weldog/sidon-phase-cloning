import numpy as np, pytest
from qhlab import designs, stabilizer, attacks, theorems, phase_hash

S1 = stabilizer.stabilizer_states(1)
S2 = stabilizer.stabilizer_states(2)

def test_theorem_a_equals_sdp_where_sw_breaks():
    for T, k, val in [(1, 3, 5 / 12), (3, 1, 5 / 6)]:
        assert abs(theorems.theorem_a(S1, T, k) - val) < 1e-9
        assert abs(attacks.counterfeit_sdp(S1, T, k) - val) < 1e-4
        assert val > designs.werner(2, T, k) + 1e-3          # S-W is false here

def test_certificate_b_is_werner_iff_t_plus_k_le_3():
    for T, k in [(1, 1), (1, 2), (2, 1)]:
        assert abs(theorems.certificate_b(S2, T, k) - designs.werner(4, T, k)) < 1e-12
    assert theorems.certificate_b(S2, 1, 3) > designs.werner(4, 1, 3) + 1e-3

def test_witting40():
    V = theorems.witting40()
    assert len(V) == 40
    assert abs(phase_hash.delta_bruteforce(V) ** 2 - 1 / 3) < 1e-9
    assert designs.is_t_design(V, 3) and not designs.is_t_design(V, 4)

def test_theorem_c_floor():
    rng = np.random.default_rng(1)
    X = np.array([[0, 1], [1, 0]]); Z = np.diag([1, -1]); I = np.eye(2)
    paulis = [np.kron(a, b) for a in (I, X, Z, X @ Z) for b in (I, X, Z, X @ Z)][1:]
    for _ in range(50):
        psi = rng.normal(size=4) + 1j * rng.normal(size=4); psi /= np.linalg.norm(psi)
        for P in paulis:
            P = P if np.allclose(P, P.conj().T) else 1j * P
            e = np.real(psi.conj() @ P @ psi)
            if abs(e) < 1 - 1e-9:
                V = (np.eye(4) + 1j * P) / np.sqrt(2)
                assert abs(psi.conj() @ V @ psi) ** 2 >= 0.5 - 1e-12

@pytest.mark.parametrize("L,K", [(6, [0, 1, 5, 17]), (7, [0, 1, 4, 9, 11])])
def test_sidon_phase_counterfeit_formula(L, K):
    d = len(K)
    assert theorems.is_sidon(K, 2 ** L)
    assert abs(attacks.counterfeit_sdp(phase_hash.codebook(L, K), 1, 1) - (2 * d - 1) / d ** 2) < 1e-4

def test_non_sidon_phase_leaks_more():
    K = [0, 1, 2, 3]
    assert not theorems.is_sidon(K, 64)
    assert attacks.counterfeit_sdp(phase_hash.codebook(6, K), 1, 1) > 7 / 16 + 0.1
