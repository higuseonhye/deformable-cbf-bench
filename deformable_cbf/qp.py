"""Equality-set solver for min 1/2 ||u - u_nom||^2 subject to A u <= b.

The decision variable is a 3-vector and the constraint count is the number
of tissues, so the active sets can be enumerated. A projected dual ascent
is kept as a check and as a fallback when M is too large to enumerate.
The dual update ascends g(lambda) = -1/2 lambda^T A A^T lambda + lambda^T (A u_nom - b).
"""

from __future__ import annotations

from itertools import combinations

import numpy as np

Array = np.ndarray


def solve_qp(u_nom: Array, A: Array, b: Array, tol: float = 1e-8) -> tuple[Array, Array]:
    """Return (u, lambda). lambda is zero when u_nom is already admissible."""
    u_nom = np.asarray(u_nom, dtype=float).reshape(-1)
    A = np.zeros((0, u_nom.shape[0])) if np.size(A) == 0 else np.asarray(A, dtype=float).reshape(-1, u_nom.shape[0])
    b = np.asarray(b, dtype=float).reshape(-1)
    count = A.shape[0]
    multipliers = np.zeros(count)
    if count == 0 or np.all(A @ u_nom <= b + tol):
        return u_nom.copy(), multipliers
    if count > 12:
        return dual_ascent(u_nom, A, b)

    for size in range(1, count + 1):
        for active in combinations(range(count), size):
            solved = _solve_active(u_nom, A, b, active, tol)
            if solved is not None:
                return solved
    return dual_ascent(u_nom, A, b)


def dual_ascent(u_nom: Array, A: Array, b: Array, iterations: int = 200) -> tuple[Array, Array]:
    """Projected gradient ascent on the dual. Step size is 1 / ||A A^T||."""
    u_nom = np.asarray(u_nom, dtype=float).reshape(-1)
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).reshape(-1)
    violation = A @ u_nom - b
    gram = A @ A.T
    lipschitz = float(np.linalg.eigvalsh(gram)[-1]) if gram.size else 1.0
    step = 1.0 / max(lipschitz, 1e-8)
    multipliers = np.zeros(A.shape[0])
    for _ in range(iterations):
        gradient = violation - gram @ multipliers
        multipliers = np.maximum(0.0, multipliers + step * gradient)
    return u_nom - A.T @ multipliers, multipliers


def _solve_active(
    u_nom: Array,
    A: Array,
    b: Array,
    active: tuple[int, ...],
    tol: float,
) -> tuple[Array, Array] | None:
    rows = list(active)
    block = A[rows]
    bound = b[rows]
    gram = block @ block.T
    rhs = block @ u_nom - bound
    multipliers_active = _solve_gram(gram, rhs)
    if multipliers_active is None:
        return None
    if np.any(multipliers_active < -1e-7):
        return None
    command = u_nom - block.T @ multipliers_active
    if np.any(A @ command - b > 1e-6):
        return None
    if np.any(np.abs(block @ command - bound) > 1e-5):
        return None
    multipliers = np.zeros(A.shape[0])
    multipliers[rows] = np.maximum(multipliers_active, 0.0)
    if np.linalg.norm(command - u_nom + A.T @ multipliers) > max(1e-5, 100 * tol):
        return None
    return command, multipliers


def _solve_gram(gram: Array, rhs: Array) -> Array | None:
    try:
        solved = np.linalg.solve(gram, rhs)
    except np.linalg.LinAlgError:
        solved, _, rank, _ = np.linalg.lstsq(gram, rhs, rcond=1e-10)
        if rank < gram.shape[0]:
            if np.linalg.norm(gram @ solved - rhs) > 1e-6 * (1.0 + np.linalg.norm(rhs)):
                return None
    else:
        if np.linalg.norm(gram @ solved - rhs) > 1e-5 * (1.0 + np.linalg.norm(rhs)):
            return None
    return solved
