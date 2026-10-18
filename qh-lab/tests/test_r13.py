import itertools, numpy as np
from qhlab import lemma_r as LR

def test_fusion_protocol_exact_for_every_xy():
    K, q = [0, 1, 2, 3], 11; d = len(K); w = np.exp(2j * np.pi / q); idx = {k: i for i, k in enumerate(K)}
    Kr = []
    for v in sorted({(b - a) % q for a in K for b in K}):
        A = np.zeros((d, d * d), complex)
        for j in range(d):
            l = idx.get((K[j] + v) % q)
            if l is not None: A[j, j * d + l] = 1
        Kr.append(A)
    assert np.allclose(sum(A.conj().T @ A for A in Kr), np.eye(d * d))
    psi = lambda x: np.array([w ** (x * k) for k in K]) / np.sqrt(d)
    E = sum(1 for a, b, c, e in itertools.product(K, repeat=4) if (a + e - b - c) % q == 0)
    for x in range(q):
        for y in range(q):
            f = sum(abs(psi((x + y) % q).conj() @ (A @ np.kron(psi(x), psi(y)))) ** 2 for A in Kr)
            assert abs(f - E / d ** 3) < 1e-12

def test_neg_table_row_d4_12_and_d2_exception():
    assert abs(LR.lemma_r([13, 20, 21, 31], 32, 1, 2) - 0.295129) < 1e-6      # Sidon-only, d=4 (1,2)
    assert abs(LR.lemma_r([0, 1, 2, 3], 1048576, 1, 2) - 0.686077) < 1e-6     # AP worst
    K, q = [0, 26], 39                                                       # d=2, Sidon, not B3
    assert LR.is_Bh(K, q, 2) and not LR.is_Bh(K, q, 3)
    assert abs(LR.lemma_r(K, q, 1, 2) - 0.75) < 1e-8                          # attains the torus value
