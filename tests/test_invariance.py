import numpy as np

from deformable_cbf.bench import hover
from deformable_cbf.filter import TVCBFFilter
from deformable_cbf.geometry import sphere_clearance_mm, sphere_state


def test_static_hover_is_entered_and_time_varying_hover_is_not():
    static = hover(predict_motion=False)
    moving = hover(predict_motion=True)
    assert static.violated
    assert static.min_clearance_mm < -0.1
    assert not moving.violated
    assert moving.min_clearance_mm >= -1e-3
    assert moving.min_clearance_mm > static.min_clearance_mm


def test_expanding_sphere_is_not_entered():
    filt = TVCBFFilter(alpha=8.0, dt=0.01, u_max=40.0)
    position = np.array([0.0, 0.0, 12.0])
    center = np.array([0.0, 0.0, 0.0])
    radius = 5.0
    margin = 1.0
    rate = 8.0  # mm/s, the safe radius grows toward the tool
    for _ in range(30):
        state = sphere_state(
            center.reshape(1, 3),
            np.zeros((1, 3)),
            np.array([radius]),
            margin,
            radius_rate=np.array([rate]),
        )
        # Nominal command dives into the sphere.
        result = filt.apply(position, np.array([0.0, 0.0, -40.0]), state, predict_motion=True)
        assert result.feasible
        position = position + 0.01 * result.velocity
        radius = radius + 0.01 * rate
        assert sphere_clearance_mm(position, center, radius, margin) >= -1e-3


def test_two_spheres_keep_a_dive_outside_both():
    filt = TVCBFFilter(alpha=8.0, dt=0.01, u_max=40.0)
    centers = np.array([[-14.0, 0.0, 0.0], [14.0, 0.0, 0.0]])
    state = sphere_state(centers, np.zeros_like(centers), np.array([5.0, 5.0]), margin=1.0)
    # 0.2 mm outside the left expanded sphere, diving toward it.
    position = np.array([-14.0, 0.0, 6.2])
    result = filt.apply(position, np.array([0.0, 0.0, -40.0]), state, predict_motion=True)
    nxt = position + 0.01 * result.velocity
    assert result.feasible
    assert result.velocity[2] > -20.0
    assert sphere_clearance_mm(nxt, centers[0], 5.0, 1.0) >= -1e-3
    assert sphere_clearance_mm(nxt, centers[1], 5.0, 1.0) >= -1e-3
