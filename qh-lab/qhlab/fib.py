"""Fibonacci anyon braid representation (benchmark generator only -- KILLED as a hash family in Round 3).
Fusion-tree basis: x_0=1, x_1=t, x_k in x_{k-1} (x) t, x_n = total charge. Labels: 0 = vacuum, 1 = tau."""
import itertools
import numpy as np

PHI = (1 + 5 ** 0.5) / 2
F = np.array([[1 / PHI, PHI ** -0.5], [PHI ** -0.5, -1 / PHI]])
R = {0: np.exp(-4j * np.pi / 5), 1: np.exp(3j * np.pi / 5)}

def basis(n, total):
    out = []
    for mid in itertools.product([0, 1], repeat=max(n - 2, 0)):
        seq = (0, 1) + mid + (total,) if n >= 2 else None
        if seq is None:
            continue
        ok = all(not (seq[i - 1] == 0 and seq[i] == 0) for i in range(1, len(seq)))
        if ok and len(seq) == n + 1:
            out.append(seq)
    return out

def generators(n, total):
    B = basis(n, total)
    idx = {b: i for i, b in enumerate(B)}
    dim = len(B)
    gens = []
    for i in range(1, n):
        G = np.zeros((dim, dim), complex)
        for b in B:
            if i == 1:
                G[idx[b], idx[b]] = R[b[2]]
                continue
            a, xi, c = b[i - 1], b[i], b[i + 1]
            if a == 1 and c == 1:
                # 2-dim space: x_i in {1, t}; channel basis via F
                for xo in (0, 1):
                    b2 = b[:i] + (xo,) + b[i + 1:]
                    if b2 in idx:
                        G[idx[b2], idx[b]] = sum(F[xo, ch] * R[ch] * F[ch, xi] for ch in (0, 1))
            else:
                ch = 0 if (a == 0 and c == 0) else 1
                G[idx[b], idx[b]] = R[ch]
        gens.append(G)
    return gens

def braid_state(word, gens, start=0):
    v = np.zeros(gens[0].shape[0], complex); v[start] = 1
    for g in word:
        i, s = abs(g) - 1, 1 if g > 0 else -1
        v = (gens[i] if s > 0 else gens[i].conj().T) @ v
    return v
