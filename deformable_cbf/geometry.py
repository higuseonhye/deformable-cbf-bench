"""Ellipsoidal barrier geometry.

The safe side of obstacle i is the exterior of a margin-expanded ellipsoid:

    h_i(x) = (x - p_i)^T Q_i (x - p_i) - 1

Q_i already contains the 1 mm margin, so h_i >= 0 is the kept set. For a
sphere of anatomical radius r and margin m, Q = I / (r + m)^2.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class EllipsoidState:
    """M obstacles. Lengths are millimetres and millimetres per second."""

    center: Array  # (M, 3)
    velocity: Array  # (M, 3)
    Q: Array  # (M, 3, 3)
    Q_dot: Array  # (M, 3, 3)

    def __post_init__(self) -> None:
        object.__setattr__(self, "center", np.asarray(self.center, dtype=float))
        object.__setattr__(self, "velocity", np.asarray(self.velocity, dtype=float))
        object.__setattr__(self, "Q", np.asarray(self.Q, dtype=float))
        object.__setattr__(self, "Q_dot", np.asarray(self.Q_dot, dtype=float))
        m = self.center.shape[0]
        if self.center.shape != (m, 3):
            raise ValueError("center must be (M, 3)")
        if self.velocity.shape != (m, 3):
            raise ValueError("velocity must be (M, 3)")
        if self.Q.shape != (m, 3, 3) or self.Q_dot.shape != (m, 3, 3):
            raise ValueError("Q and Q_dot must be (M, 3, 3)")


def Q_from_axes(radii: Array, rotation: Array) -> Array:
    """Q = R diag(1/r^2) R^T. `rotation` columns are the principal axes."""
    radii = np.asarray(radii, dtype=float)
    rotation = np.asarray(rotation, dtype=float)
    inverse = 1.0 / np.maximum(radii, 1e-6) ** 2
    return (rotation * inverse) @ rotation.T


def sphere_state(
    center: Array,
    velocity: Array,
    radius: Array,
    margin: float,
    radius_rate: Array | None = None,
) -> EllipsoidState:
    """One Q per sphere. The margin is inside Q, not added again later."""
    center = np.atleast_2d(np.asarray(center, dtype=float))
    velocity = np.atleast_2d(np.asarray(velocity, dtype=float))
    radius = np.asarray(radius, dtype=float).reshape(-1)
    count = center.shape[0]
    if radius_rate is None:
        radius_rate = np.zeros(count)
    radius_rate = np.asarray(radius_rate, dtype=float).reshape(-1)
    safe = np.maximum(radius + margin, 1e-6)
    q_scale = 1.0 / safe**2
    q_dot_scale = -2.0 * radius_rate / safe**3
    eye = np.eye(3)
    Q = q_scale[:, None, None] * eye
    Q_dot = q_dot_scale[:, None, None] * eye
    return EllipsoidState(center, velocity, Q, Q_dot)


def lie_constraints(
    position: Array,
    state: EllipsoidState,
    alpha: float,
    predict_motion: bool,
) -> tuple[Array, Array, Array]:
    """Return h (M,), A (M, 3), b (M,) for A u <= b.

    The inequality is dh/dt + alpha h >= 0 with u = dx/dt. If predict_motion
    is false, tissue velocity and Q_dot are dropped. That is the static
    barrier: it sees the current shape and ignores that the shape is moving.
    """
    diff = np.asarray(position, dtype=float).reshape(1, 3) - state.center
    q_diff = np.einsum("mij,mj->mi", state.Q, diff)
    h = np.einsum("mi,mi->m", diff, q_diff) - 1.0
    A = -2.0 * q_diff
    if predict_motion:
        velocity = state.velocity
        q_dot = state.Q_dot
    else:
        velocity = np.zeros_like(state.velocity)
        q_dot = np.zeros_like(state.Q_dot)
    transported = -2.0 * np.einsum("mi,mi->m", q_diff, velocity)
    shape_rate = np.einsum("mi,mi->m", diff, np.einsum("mij,mj->mi", q_dot, diff))
    b = alpha * h + transported + shape_rate
    return h, A, b


def predict_h(
    position: Array,
    velocity: Array,
    state: EllipsoidState,
    dt: float,
    predict_motion: bool,
) -> Array:
    """Barrier values at x + dt u under one Euler step of the tissue."""
    if predict_motion:
        center = state.center + dt * state.velocity
        Q = state.Q + dt * state.Q_dot
    else:
        center = state.center
        Q = state.Q
    Q = 0.5 * (Q + np.swapaxes(Q, -1, -2))
    diff = np.asarray(position, dtype=float).reshape(1, 3) + dt * np.asarray(velocity, dtype=float).reshape(1, 3) - center
    q_diff = np.einsum("mij,mj->mi", Q, diff)
    return np.einsum("mi,mi->m", diff, q_diff) - 1.0


def sphere_clearance_mm(position: Array, center: Array, radius: float, margin: float) -> float:
    """Euclidean gap outside the expanded sphere. Negative means the margin is gone."""
    return float(np.linalg.norm(np.asarray(position) - np.asarray(center)) - radius - margin)
