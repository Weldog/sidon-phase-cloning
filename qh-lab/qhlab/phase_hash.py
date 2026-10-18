"""QDH-1 phase hash (Ablayev-type): |psi_x> = d^{-1/2} sum_j w^{k_j x}|j>, w = exp(2 pi i / 2^L)."""
import numpy as np

def codebook(L, K):
    q = 2 ** L
    K = np.asarray(K)
    x = np.arange(q)
    return np.exp(2j * np.pi * np.outer(x, K) / q) / np.sqrt(len(K))

def delta_fft(L, K):
    """Max overlap delta = max_{a!=0} |(1/d) sum_j w^{k_j a}|, one FFT, O(q log q)."""
    q = 2 ** L
    ind = np.zeros(q)
    np.add.at(ind, np.asarray(K) % q, 1.0)
    f = np.fft.fft(ind)
    return float(np.max(np.abs(f[1:])) / len(K))

def delta_bruteforce(states):
    S = np.asarray(states)
    G = np.abs(S.conj() @ S.T)
    np.fill_diagonal(G, 0.0)
    return float(G.max())

def random_K(L, d, rng):
    return rng.choice(2 ** L, size=d, replace=False)

def welch_floor(N, d):
    """Welch bound on max |<a|b>|^2."""
    return (N - d) / (d * (N - 1))

def levenstein_floor(N, d):
    return (2 * N - d * (d + 1)) / ((d + 1) * (N - d)) if N > d * d else None
