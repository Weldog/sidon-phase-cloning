import itertools, random, collections, numpy as np
from qhlab import lemma_r as LR

def _fusion_parts(K, q, m):
    d = len(K); w = np.exp(2j * np.pi / q); idx = {k: i for i, k in enumerate(K)}
    psi = lambda x: np.array([w ** (x * k) for k in K]) / np.sqrt(d)
    br = collections.defaultdict(list)                       # delta -> list of j-tuples
    for j in itertools.product(range(d), repeat=m):
        br[tuple((K[j[i]] - K[j[0]]) % q for i in range(1, m))].append(j)
    Kraus = []
    for js in br.values():
        A = np.zeros((d, d ** m), complex)
        for j in js: A[j[0], int(np.ravel_multi_index(j, (d,) * m))] = 1
        Kraus.append(A)
    N = {dl: len(js) for dl, js in br.items()}
    r = collections.Counter((b - a) % q for a in K for b in K); Em = sum(c ** m for c in r.values())
    return d, psi, Kraus, N, Em

def test_mfold_fusion_m3_primal_and_dual():
    K, q, m = [0, 1, 3], 7, 3
    d, psi, Kraus, N, Em = _fusion_parts(K, q, m)
    assert np.allclose(sum(A.conj().T @ A for A in Kraus), np.eye(d ** m))
    X = np.zeros((d ** m * d,) * 2, complex)
    for xs in itertools.product(range(q), repeat=m):
        vin = np.array([1.0 + 0j])
        for x in xs: vin = np.kron(vin, psi(x))
        f = sum(abs(psi(sum(xs) % q).conj() @ (A @ vin)) ** 2 for A in Kraus)
        assert abs(f - Em / d ** (m + 1)) < 1e-12                      # every input, primal
        v = np.kron(vin.conj(), psi(sum(xs) % q)); X += np.outer(v, v.conj())
    X /= q ** m
    Y = np.diag([N[tuple((K[j[i]] - K[j[0]]) % q for i in range(1, m))] / d ** (m + 1)
                 for j in itertools.product(range(d), repeat=m)])
    assert np.linalg.eigvalsh(np.kron(Y, np.eye(d)) - X).min() > -1e-12   # dual feasible
    assert abs(np.trace(Y).real - Em / d ** (m + 1)) < 1e-12 and Em == 33   # = 33/81

def test_t1_exact_characterisation_random():
    def lam(K, q, k):
        r = LR.rep_counts_sparse(K, q, k + 1); return max(sum(r.get((s + a) % q, 0) for a in K) for s in range(q))
    random.seed(11)
    for _ in range(25):
        d = random.choice([3, 4]); k = random.choice([1, 2]); q = random.randint(d + 2, 30)
        K = sorted(random.sample(range(q), d)); Kt, qt = LR.torus_K(d, 1, k)
        pred = lam(K, q, k) <= lam(Kt, qt, k)
        assert pred == (abs(LR.lemma_r(K, q, 1, k) - LR.lemma_r(Kt, qt, 1, k)) < 1e-7)

def test_sidon_not_B3_attains_torus():
    K, q = [0, 1, 13], 39
    assert LR.is_Bh(K, q, 2) and not LR.is_Bh(K, q, 3)
    assert abs(LR.lemma_r(K, q, 1, 2) - 4 / 9) < 1e-8
