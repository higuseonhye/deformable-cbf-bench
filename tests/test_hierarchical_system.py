import numpy as np

from deformable_cbf.core.adaptive_phase_cbf import AdaptivePhaseCBFManager, SurgicalPhase
from deformable_cbf.core.biomechanical_predictor import (
    BiomechanicalForwardPredictor,
    VisionMeasurement,
    oscillator_step,
)
from deformable_cbf.core.box_constrained_qp import BoxConstrainedTVCBFSolver


def test_full_hierarchical_loop():
    phase_mgr = AdaptivePhaseCBFManager()
    predictor = BiomechanicalForwardPredictor(breathing_freq_hz=0.3)
    qp_solver = BoxConstrainedTVCBFSolver(u_max=0.05)

    meas = VisionMeasurement(
        timestamp=0.0,
        center=np.array([0.0, 0.0, 0.02]),
        radius=0.004,
        confidence=0.9,
    )
    predictor.update_vision(meas)
    phase_mgr.transition_phase(SurgicalPhase.HEMOSTASIS)

    t_control = 0.05
    p_pred, v_pred, uncert = predictor.predict(t_control)
    params = phase_mgr.get_interpolated_params(dt=0.01)

    ee_pos = p_pred + np.array([0.0, 0.0, 0.0055])
    u_nom = np.array([0.0, 0.0, -0.06])
    safe_radius = 0.004 + params.d_min + uncert
    diff = ee_pos - p_pred
    h = np.dot(diff, diff) - (safe_radius**2)
    A_cbf = (-2.0 * diff).reshape(1, 3)
    b_cbf = np.array([params.alpha * h - 2.0 * np.dot(diff, v_pred)])

    u_safe, intervened, _saturated = qp_solver.solve(u_nom, A_cbf, b_cbf)

    assert np.all(np.abs(u_safe) <= 0.05 + 1e-6)
    assert intervened is True
    assert float(u_safe[2]) > float(u_nom[2])


def test_dropout_grows_the_margin_and_a_dim_frame_does_not_reset_it():
    predictor = BiomechanicalForwardPredictor()
    predictor.update_vision(
        VisionMeasurement(0.0, np.zeros(3), radius=0.004, confidence=0.9)
    )
    _, _, at_100_ms = predictor.predict(0.10)
    assert at_100_ms == np.float64(0.0005) or abs(at_100_ms - 0.0005) < 1e-12

    predictor.update_vision(
        VisionMeasurement(0.10, np.zeros(3), radius=0.004, confidence=0.1)
    )
    _, _, still = predictor.predict(0.10)
    _, _, at_200_ms = predictor.predict(0.20)
    assert abs(still - 0.0005) < 1e-12
    assert abs(at_200_ms - 0.001) < 1e-12


def test_velocity_carries_the_centre_across_a_dropped_frame():
    predictor = BiomechanicalForwardPredictor()
    predictor.update_vision(VisionMeasurement(0.0, np.zeros(3), 0.004, 0.9))
    predictor.update_vision(VisionMeasurement(1.0 / 30.0, np.array([0.0, 0.0, 0.001]), 0.004, 0.9))
    held, _, _ = predictor.predict(1.0 / 30.0)
    ahead, velocity, _ = predictor.predict(1.0 / 30.0 + 0.05)
    assert ahead[2] > held[2]
    assert velocity[2] > 0.0


def test_breathing_spring_pulls_an_offset_back_when_velocity_is_zero():
    deviation, velocity = oscillator_step(
        np.array([0.0, 0.0, 0.002]),
        np.zeros(3),
        tau=0.05,
        omega=2.0 * np.pi * 0.3,
        gamma=0.95,
    )
    assert deviation[2] < 0.002
    assert velocity[2] < 0.0


def test_phase_margin_shrinks_without_a_one_tick_jump():
    manager = AdaptivePhaseCBFManager()
    manager.transition_phase(SurgicalPhase.HEMOSTASIS)
    stepped = manager.get_interpolated_params(dt=0.01)
    assert 0.0003 < stepped.d_min < 0.0010
    assert stepped.f_max == 4.0
    assert abs(stepped.d_min - (0.95 * 0.001 + 0.05 * 0.0003)) < 1e-12


def test_speed_box_clips_a_command_the_barrier_does_not_touch():
    solver = BoxConstrainedTVCBFSolver(u_max=0.05)
    command, intervened, saturated = solver.solve(
        np.array([0.0, 0.0, -0.20]),
        np.zeros((0, 3)),
        np.zeros(0),
    )
    assert np.allclose(command, [0.0, 0.0, -0.05])
    assert intervened
    assert saturated
