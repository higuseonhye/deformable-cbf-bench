"""Smooth the barrier margin when the surgical phase changes.

Lengths are metres. One 10 ms tick moves 5% of the way to the new margin
and gain, so a switch lands in about 0.2 s instead of in a single command.
The force cap changes with the phase immediately: hemostasis is allowed to
push before the position margin has finished shrinking.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class SurgicalPhase(Enum):
    APPROACH = "APPROACH"
    DISSECTION = "DISSECTION"
    HEMOSTASIS = "HEMOSTASIS"
    RETRACTION = "RETRACTION"


@dataclass(frozen=True)
class PhaseCBFParams:
    d_min: float  # metres
    alpha: float
    f_max: float  # newtons


class AdaptivePhaseCBFManager:
    def __init__(self) -> None:
        self.phase_registry = {
            SurgicalPhase.APPROACH: PhaseCBFParams(d_min=0.0015, alpha=20.0, f_max=1.0),
            SurgicalPhase.DISSECTION: PhaseCBFParams(d_min=0.0010, alpha=25.0, f_max=2.5),
            SurgicalPhase.HEMOSTASIS: PhaseCBFParams(d_min=0.0003, alpha=40.0, f_max=4.0),
            SurgicalPhase.RETRACTION: PhaseCBFParams(d_min=0.0008, alpha=15.0, f_max=3.0),
        }
        self.current_phase = SurgicalPhase.DISSECTION
        self.smooth_d_min = self.phase_registry[self.current_phase].d_min
        self.smooth_alpha = self.phase_registry[self.current_phase].alpha

    def transition_phase(self, new_phase: SurgicalPhase) -> None:
        if new_phase not in self.phase_registry:
            raise KeyError(new_phase)
        self.current_phase = new_phase

    def get_interpolated_params(self, dt: float = 0.01) -> PhaseCBFParams:
        target = self.phase_registry[self.current_phase]
        blend = min(1.0, 5.0 * dt)
        self.smooth_d_min = (1.0 - blend) * self.smooth_d_min + blend * target.d_min
        self.smooth_alpha = (1.0 - blend) * self.smooth_alpha + blend * target.alpha
        return PhaseCBFParams(d_min=self.smooth_d_min, alpha=self.smooth_alpha, f_max=target.f_max)
