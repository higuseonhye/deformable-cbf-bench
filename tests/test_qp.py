import numpy as np

from deformable_cbf.qp import dual_ascent, solve_qp


def test_unconstrained_command_is_returned():
    command, multipliers = solve_qp(np.array([1.0, 0.0, 0.0]), np.array([[1.0, 0.0, 0.0]]), np.array([2.0]))
    assert np.allclose(command, [1.0, 0.0, 0.0])
    assert np.allclose(multipliers, 0.0)


def test_single_halfspace_matches_the_closed_form():
    # u_x <= 0, nominal points along +x. Projection is the origin on x.
    A = np.array([[1.0, 0.0, 0.0]])
    b = np.array([0.0])
    u_nom = np.array([1.0, 0.4, -0.2])
    command, multipliers = solve_qp(u_nom, A, b)
    assert np.allclose(command, [0.0, 0.4, -0.2], atol=1e-8)
    assert multipliers[0] > 0.0
    assert abs(float(np.squeeze(A @ command - b))) < 1e-8


def test_kkt_on_two_active_bounds():
    A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    b = np.array([0.1, 0.2])
    u_nom = np.array([1.0, 1.0, 0.5])
    command, multipliers = solve_qp(u_nom, A, b)
    assert np.allclose(command, [0.1, 0.2, 0.5], atol=1e-7)
    stationarity = command - u_nom + A.T @ multipliers
    assert np.linalg.norm(stationarity) < 1e-7
    assert np.all(A @ command - b <= 1e-7)
    assert np.all(multipliers >= -1e-8)
    assert np.all(np.abs(multipliers * (A @ command - b)) < 1e-7)


def test_dual_ascent_matches_the_active_set():
    A = np.array(
        [
            [1.0, 0.2, 0.0],
            [-0.3, 1.0, 0.1],
            [0.0, -0.4, 1.0],
        ]
    )
    b = np.array([0.05, 0.02, 0.0])
    u_nom = np.array([0.4, -0.3, 0.5])
    exact, _ = solve_qp(u_nom, A, b)
    approx, _ = dual_ascent(u_nom, A, b, iterations=400)
    assert np.allclose(approx, exact, atol=1e-3)
    assert np.all(A @ exact - b <= 1e-6)
