"""Pad each semi-axis by a 3-D Gaussian spectral radius.

kappa = sqrt(chi^2_{3}(0.999)) = 4.033142843. Adding kappa * sigma_max to
every axis is a conservative pad around the estimated centre. It is not a
tight chance constraint on h.
"""

from __future__ import annotations

import numpy as np

from .geometry import Array, EllipsoidState, Q_from_axes

KAPPA_99_9 = 4.033142843


def inflated_radii(radii: Array, covariance: Array, margin: float, kappa: float = KAPPA_99_9) -> Array:
    """radii (M, 3) millimetres, covariance (M, 3, 3) millimetres squared."""
    radii = np.asarray(radii, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    sigma = np.sqrt(np.maximum(np.linalg.eigvalsh(covariance)[:, -1], 0.0))
    return radii + margin + kappa * sigma[:, None]


def state_from_ellipsoids(
    center: Array,
    velocity: Array,
    radii: Array,
    rotation: Array,
    margin: float,
    covariance: Array | None = None,
    radii_rate: Array | None = None,
    kappa: float = KAPPA_99_9,
) -> EllipsoidState:
    """Build Q from margin-expanded, and optionally covariance-padded, axes.

    radii_rate is d(radii)/dt before padding. The pad itself is treated as
    constant over the control step, so it does not add a Q_dot term.
    """
    center = np.atleast_2d(np.asarray(center, dtype=float))
    velocity = np.atleast_2d(np.asarray(velocity, dtype=float))
    radii = np.asarray(radii, dtype=float).reshape(center.shape[0], 3)
    rotation = np.asarray(rotation, dtype=float).reshape(center.shape[0], 3, 3)
    if covariance is None:
        safe = radii + margin
    else:
        safe = inflated_radii(radii, covariance, margin, kappa)
    if radii_rate is None:
        radii_rate = np.zeros_like(radii)
    radii_rate = np.asarray(radii_rate, dtype=float).reshape(radii.shape)
    Q = np.stack([Q_from_axes(safe[i], rotation[i]) for i in range(center.shape[0])])
    # d/dt 1/(r+pad)^2 = -2 /(r+pad)^3 * dr/dt, in the principal frame.
    Q_dot = np.zeros_like(Q)
    for i in range(center.shape[0]):
        scale = -2.0 * radii_rate[i] / np.maximum(safe[i], 1e-6) ** 3
        Q_dot[i] = (rotation[i] * scale) @ rotation[i].T
    return EllipsoidState(center, velocity, Q, Q_dot)
