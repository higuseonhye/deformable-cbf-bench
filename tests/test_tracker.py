import numpy as np

from deformable_cbf.geometry import lie_constraints
from deformable_cbf.tracker import VertexEllipsoidTracker, fit_containing_ellipsoid


def _cloud() -> np.ndarray:
    rng = np.random.default_rng(0)
    axes = np.array([6.0, 2.0, 1.5])
    return rng.normal(size=(80, 3)) * axes + np.array([10.0, -4.0, 2.0])


def test_fitted_ellipsoid_contains_every_vertex_before_the_margin_is_the_only_slack():
    points = _cloud()
    center, radii, rotation = fit_containing_ellipsoid(points, margin_mm=0.0)
    coords = (points - center) @ rotation
    quad = np.sum((coords / radii) ** 2, axis=1)
    assert np.max(quad) <= 1.0 + 1e-8
    assert np.max(quad) > 0.99


def test_translation_shows_up_as_centre_velocity():
    points = _cloud()
    tracker = VertexEllipsoidTracker(margin_mm=1.0)
    tracker.update(points, dt=0.1)
    state = tracker.update(points + np.array([1.0, 0.0, 0.0]), dt=0.1)
    assert np.allclose(state.velocity[0], [10.0, 0.0, 0.0], atol=1e-6)
    h, _, _ = lie_constraints(points[0], state, alpha=8.0, predict_motion=True)
    # A vertex of the anatomical cloud is inside the margin-expanded ellipsoid.
    assert float(h[0]) < 0.0
