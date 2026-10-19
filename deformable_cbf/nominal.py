"""Nominal velocities the barrier filters. Not a trained flow field.

The attractor is a saturated proportional command. The slide adds a tangent
when that command points into the nearest ellipsoid, so a head-on target
does not stall on the outward normal. Both pieces are then passed through
the quadratic program; neither one is a safety certificate by itself.
"""

from __future__ import annotations

import numpy as np

from .geometry import Array, EllipsoidState, lie_constraints


def attractor(position: Array, target: Array, gain: float, u_max: float) -> Array:
    delta = np.asarray(target, dtype=float).reshape(3) - np.asarray(position, dtype=float).reshape(3)
    command = gain * delta
    return _saturate(command, u_max)


def slide(
    position: Array,
    command: Array,
    state: EllipsoidState,
    target: Array,
    u_max: float,
    alpha: float = 8.0,
    h_trigger: float = 1.0,
) -> Array:
    """Add a tangent toward the target when the command enters a nearby zone."""
    position = np.asarray(position, dtype=float).reshape(3)
    command = np.asarray(command, dtype=float).reshape(3)
    h, _, _ = lie_constraints(position, state, alpha, predict_motion=True)
    index = int(np.argmin(h))
    if h[index] > h_trigger:
        return command
    diff = position - state.center[index]
    outward = state.Q[index] @ diff
    norm = float(np.linalg.norm(outward))
    if norm < 1e-9:
        return command
    normal = outward / norm
    if float(normal @ command) >= 0.0:
        return command
    to_target = np.asarray(target, dtype=float).reshape(3) - position
    tangent = to_target - normal * float(normal @ to_target)
    tangent_norm = float(np.linalg.norm(tangent))
    if tangent_norm < 1e-6:
        helper = np.array([1.0, 0.0, 0.0]) if abs(normal[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        tangent = np.cross(normal, helper)
        tangent_norm = float(np.linalg.norm(tangent))
    tangent = tangent / tangent_norm
    return _saturate(command + u_max * tangent, u_max)


def replan(position: Array, state: EllipsoidState, u_max: float) -> Array:
    """Step away from the nearest centre, with a fixed lateral bias."""
    position = np.asarray(position, dtype=float).reshape(3)
    offset = position.reshape(1, 3) - state.center
    index = int(np.argmin(np.linalg.norm(offset, axis=1)))
    outward = offset[index]
    norm = float(np.linalg.norm(outward))
    normal = np.array([0.0, 0.0, 1.0]) if norm < 1e-9 else outward / norm
    lateral = np.cross(normal, np.array([0.0, 0.0, 1.0]))
    if float(np.linalg.norm(lateral)) < 1e-6:
        lateral = np.array([1.0, 0.0, 0.0])
    lateral = lateral / np.linalg.norm(lateral)
    direction = 0.6 * normal + 0.8 * lateral + np.array([0.0, 0.0, 0.35])
    return _saturate(direction, u_max)


def _saturate(command: Array, u_max: float) -> Array:
    command = np.asarray(command, dtype=float).reshape(3)
    speed = float(np.linalg.norm(command))
    if speed > u_max > 0.0:
        return command * (u_max / speed)
    return command
