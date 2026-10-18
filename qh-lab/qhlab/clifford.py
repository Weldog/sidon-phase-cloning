"""Clifford group (mod phase) for n=1,2 and Clifford-orbit codebooks.
Clifford group is a unitary 3-design => ANY fiducial orbit is a state 3-design (eta_{1,1}=0).
Collision bias delta of an orbit needs only |C| fiducial overlaps: <g f|h f> = <f|g^dag h|f>."""
import numpy as np

def _canon_u(U):
    flat = U.ravel()
    i = np.flatnonzero(np.abs(flat) > 1e-9)[0]
    U = U * (abs(flat[i]) / flat[i])
    return tuple(np.round(U.real, 8).ravel()) + tuple(np.round(U.imag, 8).ravel()), U

def clifford_group(n):
    from .stabilizer import stabilizer_states  # noqa (same gate set)
    dim = 2 ** n
    H1 = np.array([[1, 1], [1, -1]]) / np.sqrt(2); S1 = np.diag([1, 1j])
    def one(g, q):
        out = np.array([[1.0 + 0j]])
        for j in range(n):
            out = np.kron(out, g if j == q else np.eye(2))
        return out
    gens = [one(H1, q) for q in range(n)] + [one(S1, q) for q in range(n)]
    if n == 2:
        gens.append(np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]], complex))
    key, I = _canon_u(np.eye(dim, dtype=complex))
    seen = {key: I}; frontier = [I]
    while frontier:
        nxt = []
        for U in frontier:
            for g in gens:
                k2, V = _canon_u(g @ U)
                if k2 not in seen:
                    seen[k2] = V; nxt.append(V)
        frontier = nxt
    return np.array(list(seen.values()))

def orbit_delta(group, f):
    ov = np.abs(np.einsum('i,gij,j->g', f.conj(), group, f))
    return float(ov[ov < 1 - 1e-9].max())

def orbit_states(group, f):
    S = group @ f
    keys, out = set(), []
    for s in S:
        i = np.flatnonzero(np.abs(s) > 1e-9)[0]
        s = s * (abs(s[i]) / s[i])
        k = tuple(np.round(np.concatenate([s.real, s.imag]), 8))
        if k not in keys:
            keys.add(k); out.append(s)
    return np.array(out)

def best_fiducial(group, trials, rng):
    dim = group.shape[1]
    best, bf = 2.0, None
    for _ in range(trials):
        f = rng.normal(size=dim) + 1j * rng.normal(size=dim); f /= np.linalg.norm(f)
        dlt = orbit_delta(group, f)
        if dlt < best:
            best, bf = dlt, f
    # local refinement
    step = 0.05
    for _ in range(4000):
        g = bf + step * (rng.normal(size=dim) + 1j * rng.normal(size=dim)); g /= np.linalg.norm(g)
        dlt = orbit_delta(group, g)
        if dlt < best:
            best, bf = dlt, g
        else:
            step *= 0.999
    return best, bf
