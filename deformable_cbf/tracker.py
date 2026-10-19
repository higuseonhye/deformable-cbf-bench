"""Fit one containing ellipsoid to a vertex cluster.

The axes are the principal directions of the point cloud. Radii are scaled
so every vertex satisfies (x-p)^T Q (x-p) <= 1 before the safety margin is
added. Rates are a backward difference of the fitted centre and of Q.
Isaac is not required: pass an (N, 3) array from any mesh buffer.
"""

from __future__ import annotations

import numpy as np

from .geometry import Array, EllipsoidState, Q_from_axes


class VertexEllipsoidTracker:
    def __init__(self, margin_mm: float = 1.0):
        self.margin_mm = float(margin_mm)
        self._previous_center: Array | None = None
        self._previous_Q: Array | None = None

    def update(self, vertices: Array, dt: float) -> EllipsoidState:
        center, radii, rotation = fit_containing_ellipsoid(vertices, self.margin_mm)
        Q = Q_from_axes(radii, rotation)[None, ...]
        if self._previous_center is None or dt <= 0.0:
            velocity = np.zeros((1, 3))
            Q_dot = np.zeros((1, 3, 3))
        else:
            velocity = (center - self._previous_center).reshape(1, 3) / dt
            Q_dot = (Q - self._previous_Q) / dt
        self._previous_center = center.copy()
        self._previous_Q = Q.copy()
        return EllipsoidState(center.reshape(1, 3), velocity, Q, Q_dot)


def fit_containing_ellipsoid(vertices: Array, margin_mm: float) -> tuple[Array, Array, Array]:
    """Return center (3,), radii including margin (3,), and rotation (3, 3)."""
    points = _as_points(vertices)
    center = points.mean(axis=0)
    centered = points - center
    _, _, vh = np.linalg.svd(centered, full_matrices=True)
    rotation = vh[:3].T
    coords = centered @ rotation
    radii = np.maximum(np.max(np.abs(coords), axis=0), 1e-3)
    scale = float(np.sqrt(np.max(np.sum((coords / radii) ** 2, axis=1))))
    radii = radii * max(scale, 1.0) + margin_mm
    return center, radii, rotation


def _as_points(vertices: Array) -> Array:
    if hasattr(vertices, "detach"):
        vertices = vertices.detach().cpu().numpy()
    points = np.asarray(vertices, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] < 4:
        raise ValueError("vertices must have shape (N, 3) with N >= 4")
    return points
