import numpy as np

from deformable_cbf.filter import TVCBFFilter
from deformable_cbf.force import force_row
from deformable_cbf.geometry import sphere_state
from deformable_cbf.robust import inflated_radii, state_from_ellipsoids


def test_force_row_cuts_the_inward_speed():
    force = np.array([0.0, 0.0, 2.99])
    _, A, b = force_row(force, f_max=3.0, stiffness=0.5, alpha=30.0)
    filt = TVCBFFilter(alpha=8.0, dt=0.01, u_max=40.0)
    state = sphere_state(np.array([[0.0, 0.0, -50.0]]), np.zeros((1, 3)), np.array([2.0]), margin=1.0)
    result = filt.apply(
        np.zeros(3),
        np.array([0.0, 0.0, 40.0]),
        state,
        predict_motion=True,
        extra_A=A,
        extra_b=b,
    )
    assert float(result.velocity[2]) < 1.0
    assert float(np.squeeze(A @ result.velocity - b)) <= 1e-6


def test_covariance_pad_is_larger_than_the_geometric_margin():
    radii = np.array([[4.0, 5.0, 3.0]])
    small = np.repeat(np.eye(3)[None, ...], 1, axis=0) * (0.1**2)
    large = np.repeat(np.eye(3)[None, ...], 1, axis=0) * (0.8**2)
    padded_small = inflated_radii(radii, small, margin=1.0)
    padded_large = inflated_radii(radii, large, margin=1.0)
    assert np.all(padded_large > padded_small)
    assert np.all(padded_small > radii + 1.0)


def test_uncertain_centre_pushes_when_the_nominal_sphere_would_allow_a_hold():
    radii = np.array([[4.0, 4.0, 4.0]])
    rotation = np.repeat(np.eye(3)[None, ...], 1, axis=0)
    covariance = np.repeat(np.eye(3)[None, ...], 1, axis=0) * (0.8**2)
    nominal = state_from_ellipsoids(np.zeros((1, 3)), np.zeros((1, 3)), radii, rotation, margin=1.0)
    robust = state_from_ellipsoids(
        np.zeros((1, 3)),
        np.zeros((1, 3)),
        radii,
        rotation,
        margin=1.0,
        covariance=covariance,
    )
    # 7 mm from the centre: outside the 5 mm expanded sphere, inside the padded one.
    position = np.array([7.0, 0.0, 0.0])
    filt = TVCBFFilter(alpha=8.0, dt=0.01, u_max=40.0)
    held = filt.apply(position, np.zeros(3), nominal, predict_motion=False)
    pushed = filt.apply(position, np.zeros(3), robust, predict_motion=False)
    assert float(np.linalg.norm(held.velocity)) < 1e-6
    assert float(pushed.velocity[0]) > 1.0
