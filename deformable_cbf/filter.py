"""Continuous barrier projection, then a one-step repair.

The half-space A u <= b is the instantaneous condition dh/dt + alpha h >= 0.
A finite control period can still step into the zone, and slowing that
command down makes the miss worse when the tissue is the thing that moves.
The repair adds speed along the outward normal of the most violated
predicted barrier, up to u_max. If the tissue outruns u_max, the margin
cannot be kept; the result records that.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import EllipsoidState, lie_constraints, predict_h
from .qp import solve_qp

Array = np.ndarray


@dataclass(frozen=True)
class FilterResult:
    velocity: Array
    intervened: bool
    h: Array
    h_next: Array
    feasible: bool


class TVCBFFilter:
    def __init__(self, alpha: float = 8.0, dt: float = 0.01, u_max: float = 40.0):
        self.alpha = float(alpha)
        self.dt = float(dt)
        self.u_max = float(u_max)

    def apply(
        self,
        position: Array,
        u_nom: Array,
        state: EllipsoidState,
        predict_motion: bool = True,
        extra_A: Array | None = None,
        extra_b: Array | None = None,
    ) -> FilterResult:
        position = np.asarray(position, dtype=float).reshape(3)
        u_nom = np.asarray(u_nom, dtype=float).reshape(3)
        h, A, b = lie_constraints(position, state, self.alpha, predict_motion)
        if extra_A is not None and np.size(extra_A):
            A = np.vstack([A, np.asarray(extra_A, dtype=float).reshape(-1, 3)])
            b = np.concatenate([b, np.asarray(extra_b, dtype=float).reshape(-1)])
        command, _ = solve_qp(u_nom, A, b)
        command = _cap(command, self.u_max)
        command, feasible = self._repair(position, command, state, predict_motion)
        h_next = predict_h(position, command, state, self.dt, predict_motion)
        intervened = float(np.linalg.norm(command - u_nom)) > 1e-6
        return FilterResult(command, intervened, h, h_next, feasible)

    def _repair(
        self,
        position: Array,
        command: Array,
        state: EllipsoidState,
        predict_motion: bool,
    ) -> tuple[Array, bool]:
        command = np.array(command, dtype=float)
        for _ in range(32):
            h_next = predict_h(position, command, state, self.dt, predict_motion)
            if np.all(h_next >= -1e-8):
                return _cap(command, self.u_max), True
            index = int(np.argmin(h_next))
            direction = _outward(position, command, state, index, self.dt, predict_motion)
            command = command + direction * max(self.u_max * 0.2, 1.0)
            command = _cap(command, self.u_max)
        h_next = predict_h(position, command, state, self.dt, predict_motion)
        return command, bool(np.all(h_next >= -1e-6))


def _outward(
    position: Array,
    command: Array,
    state: EllipsoidState,
    index: int,
    dt: float,
    predict_motion: bool,
) -> Array:
    if predict_motion:
        center = state.center[index] + dt * state.velocity[index]
        Q = state.Q[index] + dt * state.Q_dot[index]
    else:
        center = state.center[index]
        Q = state.Q[index]
    Q = 0.5 * (Q + Q.T)
    diff = position + dt * command - center
    direction = Q @ diff
    norm = float(np.linalg.norm(direction))
    if norm < 1e-9:
        return np.array([0.0, 0.0, 1.0])
    return direction / norm


def _cap(command: Array, u_max: float) -> Array:
    speed = float(np.linalg.norm(command))
    if speed > u_max > 0.0:
        return command * (u_max / speed)
    return command
