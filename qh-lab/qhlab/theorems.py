"""Round-4 results as code. Theorem A (exact Gram formula for unitary-3-group orbits, T<=3),
Theorem B certificate C_{T,k} = D_T * lambda_max(E psi^{(x)(T+k)}), Theorem C collision floor."""
import numpy as np
from .designs import D, kron_pow

def theorem_a(states, T, k):
    """(D_T/N) * lambda_max(G), G_xy = |g|^{2T} g^k. Upper bound for ANY codebook; exact for orbits of a group
    acting irreducibly on Sym_T (unitary 3-groups, T<=3)."""
    S = np.asarray(states); N, d = S.shape
    g = S.conj() @ S.T
    G = np.abs(g) ** (2 * T) * g ** k
    lam = np.linalg.eigvalsh((G + G.conj().T) / 2)[-1]
    return float(D(d, T) / N * lam)

def certificate_b(states, T, k):
    """C_{T,k} = D_T * lambda_max(M_{T+k}); equals Werner W iff codebook is a (T+k)-design. Gram trick on the
    (T+k)-fold powers keeps it cheap: nonzero spec(M_m) = spec(Gram_m)/N with Gram_m = g^m."""
    S = np.asarray(states); N, d = S.shape
    g = S.conj() @ S.T
    Gm = g ** (T + k)
    lam = np.linalg.eigvalsh((Gm + Gm.conj().T) / 2)[-1]
    return float(D(d, T) * lam / N)

def witting40():
    """Search sign conventions for 40 rays in C^4: 4 basis vectors + 36 of type (0,1,s1 w^a, s2 w^b)/sqrt3 with
    the zero cycled through 4 positions. Return the first configuration whose overlaps are only {0, 1/3}."""
    w = np.exp(2j * np.pi / 3)
    import itertools
    for signs in itertools.product([1, -1], repeat=8):
        V = [np.eye(4)[i] for i in range(4)]
        for p in range(4):
            s1, s2 = signs[2 * p], signs[2 * p + 1]
            others = [i for i in range(4) if i != p]
            for a in range(3):
                for b in range(3):
                    v = np.zeros(4, complex)
                    v[others[0]] = 1; v[others[1]] = s1 * w ** a; v[others[2]] = s2 * w ** b
                    V.append(v / np.sqrt(3))
        V = np.array(V)
        o = np.abs(V.conj() @ V.T) ** 2
        np.fill_diagonal(o, 0)
        vals = np.unique(np.round(o, 9))
        if set(vals.tolist()) <= {0.0, round(1 / 3, 9)}:
            return V
    return None

def is_sidon(K, q):
    K = [k % q for k in K]
    sums = {}
    for i in range(len(K)):
        for j in range(i, len(K)):
            s = (K[i] + K[j]) % q
            if s in sums and sums[s] != (i, j):
                return False
            sums[s] = (i, j)
    return True
