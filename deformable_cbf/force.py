"""Linear force row stacked with the geometric barrier.

The model is F = K * penetration, with K in N/mm and the command in mm/s,
so dF/dt = K u when the command is along the applied force. F is the force
the instrument applies to the tissue, in the same frame as u.

    h_F = F_max^2 - ||F||^2
    A_F = 2 K F^T
    b_F = alpha_F * h_F

A command along F increases the load. This is not a measured organ law.
"""

from __future__ import annotations

import numpy as np

from .geometry import Array


def force_row(force: Array, f_max: float, stiffness: float, alpha: float) -> tuple[float, Array, Array]:
    force = np.asarray(force, dtype=float).reshape(3)
    h = float(f_max**2 - force @ force)
    A = (2.0 * stiffness * force).reshape(1, 3)
    b = np.array([alpha * h])
    return h, A, b
