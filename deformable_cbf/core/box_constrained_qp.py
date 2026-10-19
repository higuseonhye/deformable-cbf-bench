"""Project a nominal velocity onto the barrier and a speed box.

The program is

    min  1/2 ||u - u_nom||^2
    s.t. A u <= b
         |u_i| <= u_max

Units are metres per second, matching the rest of deformable_cbf.core.
If the half-spaces and the box have no point in common, the motor box wins:
the returned command is the barrier solution clipped to u_max. The barrier
is then short of the speed it asked for, and `saturated` is true.
"""

from __future__ import annotations

import numpy as np

from deformable_cbf.qp import solve_qp

Array = np.ndarray


class BoxConstrainedTVCBFSolver:
    def __init__(self, u_max: float = 0.05):
        if u_max <= 0.0:
            raise ValueError("u_max must be positive")
        self.u_max = float(u_max)

    def solve(self, u_nom: Array, A_cbf: Array, b_cbf: Array) -> tuple[Array, bool, bool]:
        u_nom = np.asarray(u_nom, dtype=float).reshape(3)
        A_cbf = _rows(A_cbf, 3)
        b_cbf = np.asarray(b_cbf, dtype=float).reshape(-1)
        if A_cbf.shape[0] != b_cbf.shape[0]:
            raise ValueError("A_cbf and b_cbf have different row counts")

        barrier_u, _ = solve_qp(u_nom, A_cbf, b_cbf)
        if _inside_box(barrier_u, self.u_max) and _respects(A_cbf, barrier_u, b_cbf):
            command = barrier_u
        else:
            box_A = np.vstack([np.eye(3), -np.eye(3)])
            box_b = np.full(6, self.u_max)
            A_all = np.vstack([A_cbf, box_A]) if A_cbf.size else box_A
            b_all = np.concatenate([b_cbf, box_b]) if b_cbf.size else box_b
            joint, _ = solve_qp(u_nom, A_all, b_all)
            if _inside_box(joint, self.u_max) and _respects(A_cbf, joint, b_cbf):
                command = joint
            else:
                command = np.clip(barrier_u, -self.u_max, self.u_max)

        intervened = bool(np.linalg.norm(command - u_nom) > 1e-8)
        saturated = bool(np.any(np.abs(np.abs(command) - self.u_max) <= 1e-8))
        return command, intervened, saturated


def _rows(matrix: Array, width: int) -> Array:
    if np.size(matrix) == 0:
        return np.zeros((0, width))
    return np.asarray(matrix, dtype=float).reshape(-1, width)


def _inside_box(command: Array, u_max: float) -> bool:
    return bool(np.all(np.abs(command) <= u_max + 1e-8))


def _respects(A: Array, command: Array, b: Array) -> bool:
    if A.size == 0:
        return True
    return bool(np.all(A @ command - b <= 1e-6))
