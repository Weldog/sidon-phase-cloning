"""n-qubit stabilizer states via Clifford-orbit BFS (exact complex projective 3-design)."""
import numpy as np

def _canon(v):
    i = np.flatnonzero(np.abs(v) > 1e-9)[0]
    v = v * (abs(v[i]) / v[i])
    return tuple(np.round(v.real, 9)) + tuple(np.round(v.imag, 9)), v

def stabilizer_states(n):
    dim = 2 ** n
    H1 = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    S1 = np.diag([1, 1j])
    def one(g, q):
        ops = [np.eye(2)] * n
        ops[q] = g
        out = np.array([[1.0 + 0j]])
        for o in ops:
            out = np.kron(out, o)
        return out
    gates = [one(H1, q) for q in range(n)] + [one(S1, q) for q in range(n)]
    for c in range(n):
        for t in range(n):
            if c != t:
                U = np.zeros((dim, dim))
                for b in range(dim):
                    bits = [(b >> (n - 1 - i)) & 1 for i in range(n)]
                    if bits[c]:
                        bits[t] ^= 1
                    U[sum(bit << (n - 1 - i) for i, bit in enumerate(bits)), b] = 1
                gates.append(U)
    start = np.zeros(dim, complex); start[0] = 1
    key, v = _canon(start)
    seen = {key: v}
    frontier = [v]
    while frontier:
        nxt = []
        for s in frontier:
            for g in gates:
                k2, w = _canon(g @ s)
                if k2 not in seen:
                    seen[k2] = w
                    nxt.append(w)
        frontier = nxt
    return np.array(list(seen.values()))
