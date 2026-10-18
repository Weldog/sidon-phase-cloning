"""Design-theory metrics: D_T, Werner cloning value, Haar moments, design defect eta, frame potential."""
import itertools, math
import numpy as np

def D(d, T):
    """dim Sym^T(C^d) = C(d+T-1, T)."""
    return math.comb(d + T - 1, T)

def werner(d, T, k):
    """Optimal universal T->T+k cloning global fidelity (Werner 1998) = D_T / D_{T+k}."""
    return D(d, T) / D(d, T + k)

def _perm_op(d, m, perm):
    dim = d ** m
    P = np.zeros((dim, dim))
    for idx in itertools.product(range(d), repeat=m):
        src = np.ravel_multi_index(idx, (d,) * m)
        dst = np.ravel_multi_index(tuple(idx[p] for p in perm), (d,) * m)
        P[dst, src] = 1.0
    return P

def sym_projector(d, m):
    perms = list(itertools.permutations(range(m)))
    return sum(_perm_op(d, m, p) for p in perms) / len(perms)

def partial_transpose(op, d, m, sites):
    t = op.reshape((d,) * (2 * m))
    axes = list(range(2 * m))
    for s in sites:
        axes[s], axes[m + s] = axes[m + s], axes[s]
    return t.transpose(axes).reshape(d ** m, d ** m)

def kron_pow(v, p):
    out = np.array([1.0 + 0j])
    for _ in range(p):
        out = np.kron(out, v)
    return out

def moment_operator(states, T, k):
    """X = (1/N) sum_x |phi_x><phi_x|, phi_x = conj(psi_x)^{(x)T} (x) psi_x^{(x)T+k}."""
    states = np.asarray(states)
    vecs = [np.kron(kron_pow(s.conj(), T), kron_pow(s, T + k)) for s in states]
    V = np.array(vecs)
    return V.T @ V.conj() / len(V)

def haar_moment_operator(d, T, k):
    m = 2 * T + k
    return partial_transpose(sym_projector(d, m) / D(d, m), d, m, list(range(T)))

def design_defect(states, T, k):
    """eta_{T,k} = || X_codebook - X_Haar ||_op  (Round-3 definition)."""
    d = np.asarray(states).shape[1]
    diff = moment_operator(states, T, k) - haar_moment_operator(d, T, k)
    return float(np.linalg.norm(diff, 2))

def frame_potential(states, t):
    S = np.asarray(states)
    G = np.abs(S.conj() @ S.T) ** (2 * t)
    return float(G.mean())

def is_t_design(states, t, tol=1e-9):
    d = np.asarray(states).shape[1]
    return abs(frame_potential(states, t) - 1 / D(d, t)) < tol
