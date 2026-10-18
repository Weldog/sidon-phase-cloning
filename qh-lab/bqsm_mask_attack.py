"""QDH-BQSM stress test (Derived). Broadcast masked copies psi_{x+r_i} (fresh uniform shift r_i per copy, revealed after t_end).
Each masked copy alone is exactly I/d, so a memory bound looks like it should cap the attacker's copy count.
Attack with ZERO quantum memory: measure every copy at once with the shift-covariant POVM M_y=(d/q)|psi_y><psi_y|, keep
the classical outcome y_i; after the reveal, y_i - r_i is a sample centred on x. Maximum-likelihood decode."""
import json, numpy as np
rng = np.random.default_rng(20261016)

def run(K, q, T, trials=200):
    d = len(K)
    # p(y | z) depends on (y - z) only: p(u) = |sum_j w^{-u k_j}|^2 / (q d)
    ind = np.zeros(q); np.add.at(ind, np.mod(K, q), 1)
    pu = np.abs(np.fft.fft(ind)) ** 2 / (q * d)
    assert abs(pu.sum() - 1) < 1e-9          # POVM is complete (K distinct mod q)
    logp = np.log(np.maximum(pu, 1e-300))
    ok = 0
    for _ in range(trials):
        x = rng.integers(q)
        r = rng.integers(q, size=T)                       # masks, revealed after t_end
        u = rng.choice(q, size=T, p=pu)                   # outcome offsets
        y = (x + r + u) % q                               # immediate measurement outcomes (classical, stored)
        s = (y - r) % q                                   # unmask after reveal
        ll = np.array([logp[(s - c) % q].sum() for c in range(q)])
        ok += int(np.argmax(ll) == x)
    return ok / trials

out = []
for (d, L) in [(8, 10), (16, 12)]:
    q = 2 ** L
    K = np.sort(rng.choice(q, size=d, replace=False))
    for T in [2, 4, 8, 16, 32, 64]:
        out.append(dict(d=d, L=L, T=T, quantum_memory_qubits=0, success=run(K, q, T),
                        holevo_min_copies=L / np.log2(d)))
        print(out[-1], flush=True)
json.dump(out, open("bqsm_mask_attack.json", "w"), indent=1, default=float)
