"""d=2 unit test only: closest pair of a qubit codebook is a spherical-Delaunay (convex-hull) edge."""
import numpy as np
from scipy.spatial import ConvexHull

def bloch(s):
    a, b = s
    return np.array([2 * (a.conjugate() * b).real, 2 * (a.conjugate() * b).imag, abs(a) ** 2 - abs(b) ** 2])

def closest_pair_is_delaunay_edge(states):
    P = np.array([bloch(s) for s in states])
    hull = ConvexHull(P)
    edges = set()
    for simp in hull.simplices:
        for u in range(3):
            for v in range(u + 1, 3):
                edges.add(tuple(sorted((simp[u], simp[v]))))
    G = P @ P.T
    np.fill_diagonal(G, -2)
    i, j = np.unravel_index(np.argmax(G), G.shape)
    return tuple(sorted((i, j))) in edges
