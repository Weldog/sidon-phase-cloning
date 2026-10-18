"""Attack harness: optimal preimage (SDP), pretty-good measurement (Gram), optimal counterfeiting (Choi SDP)."""
import numpy as np
import cvxpy as cp
from .designs import kron_pow

SOLVER = dict(solver="SCS", eps=1e-8, max_iters=200000)

def preimage_sdp(states, T):
    """Max (1/N) sum_x Tr[M_x psi_x^{(x)T}] over POVMs. Exact optimum (small N, d^T)."""
    rhos = [np.outer(kron_pow(s, T), kron_pow(s, T).conj()) for s in states]
    D = rhos[0].shape[0]
    M = [cp.Variable((D, D), hermitian=True) for _ in rhos]
    cons = [m >> 0 for m in M] + [sum(M) == np.eye(D)]
    obj = cp.Maximize(cp.real(sum(cp.trace(M[i] @ rhos[i]) for i in range(len(rhos)))) / len(rhos))
    p = cp.Problem(obj, cons); p.solve(**SOLVER)
    return float(p.value)

def pgm_success(states, T):
    """Square-root measurement success for equiprobable pure states via Gram matrix: (1/N) sum_x (G^{1/2})_xx^2."""
    S = np.asarray(states)
    G = (S.conj() @ S.T) ** T
    w, V = np.linalg.eigh((G + G.conj().T) / 2)
    w = np.clip(w, 0, None)
    R = (V * np.sqrt(w)) @ V.conj().T
    return float(np.mean(np.abs(np.diag(R)) ** 2))

def counterfeit_sdp(states, T, k):
    """Optimal T -> T+k counterfeiting: max Tr[J X], J Choi of CPTP map (Tr_out J = I)."""
    S = np.asarray(states); d = S.shape[1]
    din, dout = d ** T, d ** (T + k)
    X = sum(np.kron(np.outer(kron_pow(s, T), kron_pow(s, T).conj()).conj(),
                    np.outer(kron_pow(s, T + k), kron_pow(s, T + k).conj())) for s in S) / len(S)
    J = cp.Variable((din * dout, din * dout), hermitian=True)
    cons = [J >> 0, cp.partial_trace(J, [din, dout], axis=1) == np.eye(din)]
    p = cp.Problem(cp.Maximize(cp.real(cp.trace(J @ X))), cons); p.solve(**SOLVER)
    return float(p.value)

def binding_bound(delta, k):
    return (1 + delta ** k) / 2
