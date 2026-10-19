"""Trajectories behind the board. Same integrator as the benchmark."""

from __future__ import annotations

import numpy as np

from .bench import ALPHA, DT, MARGIN_MM, RADIUS_MM, U_MAX, _nominal
from .filter import TVCBFFilter
from .geometry import sphere_clearance_mm, sphere_state

Array = np.ndarray


def hover_trace(predict_motion: bool, steps: int = 40, tissue_mm_s: float = 25.0) -> dict[str, Array]:
    """Tool z, expanded surface, and clearance. The command is zero."""
    filt = TVCBFFilter(alpha=ALPHA, dt=DT, u_max=U_MAX)
    position = np.array([0.0, 0.0, 8.0])
    center = np.array([0.0, 0.0, 0.0])
    velocity = np.array([0.0, 0.0, tissue_mm_s])
    times = [0.0]
    tool_z = [float(position[2])]
    surface_z = [float(center[2] + RADIUS_MM + MARGIN_MM)]
    clearance = [sphere_clearance_mm(position, center, RADIUS_MM, MARGIN_MM)]
    for step in range(steps):
        state = sphere_state(center.reshape(1, 3), velocity.reshape(1, 3), np.array([RADIUS_MM]), MARGIN_MM)
        result = filt.apply(position, np.zeros(3), state, predict_motion=predict_motion)
        position = position + DT * result.velocity
        center = center + DT * velocity
        times.append((step + 1) * DT)
        tool_z.append(float(position[2]))
        surface_z.append(float(center[2] + RADIUS_MM + MARGIN_MM))
        clearance.append(sphere_clearance_mm(position, center, RADIUS_MM, MARGIN_MM))
    return {
        "t": np.asarray(times),
        "tool_z": np.asarray(tool_z),
        "surface_z": np.asarray(surface_z),
        "clearance": np.asarray(clearance),
    }


def crossing_trace(strategy: str, phase: float = 0.0, steps: int = 300) -> dict[str, Array]:
    """Tool path in the y-z plane and clearance after each tissue step."""
    filt = TVCBFFilter(alpha=ALPHA, dt=DT, u_max=U_MAX)
    position = np.array([0.0, 18.0, 8.0])
    target = np.array([0.0, -18.0, 8.0])
    times = [0.0]
    tool = [position.copy()]
    sphere = [np.array([0.0, 0.0, 3.0 * np.sin(phase)])]
    clearance = [sphere_clearance_mm(position, sphere[0], RADIUS_MM, MARGIN_MM)]
    for step in range(steps):
        t = step * DT
        z = 3.0 * np.sin(2.0 * np.pi * 0.3 * t + phase)
        vz = 3.0 * (2.0 * np.pi * 0.3) * np.cos(2.0 * np.pi * 0.3 * t + phase)
        center = np.array([0.0, 0.0, z])
        state = sphere_state(center.reshape(1, 3), np.array([[0.0, 0.0, vz]]), np.array([RADIUS_MM]), MARGIN_MM)
        occluded = 80 <= step < 130
        command = _nominal(strategy, position, target, state, occluded)
        if strategy == "attractor":
            applied = command
        else:
            applied = filt.apply(position, command, state, predict_motion=strategy != "static").velocity
        position = position + DT * applied
        t_next = (step + 1) * DT
        z_next = 3.0 * np.sin(2.0 * np.pi * 0.3 * t_next + phase)
        center_next = np.array([0.0, 0.0, z_next])
        times.append(t_next)
        tool.append(position.copy())
        sphere.append(center_next.copy())
        clearance.append(sphere_clearance_mm(position, center_next, RADIUS_MM, MARGIN_MM))
    return {
        "t": np.asarray(times),
        "tool": np.vstack(tool),
        "sphere": np.vstack(sphere),
        "clearance": np.asarray(clearance),
        "target": target,
    }


def phase_margin_mm(steps: int = 60, dt: float = 0.01) -> tuple[Array, Array]:
    """Margin after a switch from dissection (1.0 mm) to hemostasis (0.3 mm)."""
    from .core.adaptive_phase_cbf import AdaptivePhaseCBFManager, SurgicalPhase

    manager = AdaptivePhaseCBFManager()
    manager.transition_phase(SurgicalPhase.HEMOSTASIS)
    times = np.arange(1, steps + 1) * dt
    margin = np.empty(steps)
    for index in range(steps):
        margin[index] = manager.get_interpolated_params(dt).d_min * 1000.0
    return times, margin


def dropout_margin_mm(tau_s: Array | None = None) -> tuple[Array, Array]:
    """Extra margin, in millimetres, after a vision gap. Caps at 200 ms."""
    if tau_s is None:
        tau_s = np.linspace(0.0, 0.20, 41)
    tau = np.clip(np.asarray(tau_s, dtype=float), 0.0, 0.20)
    return tau, 0.5 * (tau / 0.10)
