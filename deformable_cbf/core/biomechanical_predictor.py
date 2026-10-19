"""Bridge a 30 Hz vision sample to a 100 Hz tissue state.

Between accepted frames the centre obeys

    ÿ + γ ẏ + ω_b² y = 0

about a slow baseline. y is the deviation of the filtered centre from that
baseline, ω_b is the breathing frequency, and γ is the velocity decay in
1/s. A 10–100 ms gap is short next to a 0.3 Hz breath, so the step is close
to constant velocity; the spring term is what pulls a held offset back.

Frames with confidence at or below 0.4 do not reset the clock. Uncertainty
then grows by 0.5 mm per 100 ms of dropout, and extrapolation stops at 200 ms.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class VisionMeasurement:
    timestamp: float
    center: Array  # (3,) metres
    radius: float  # metres
    confidence: float


class BiomechanicalForwardPredictor:
    def __init__(self, breathing_freq_hz: float = 0.3, damping: float = 0.95):
        if breathing_freq_hz <= 0.0:
            raise ValueError("breathing_freq_hz must be positive")
        self.omega_b = 2.0 * np.pi * float(breathing_freq_hz)
        self.damping = float(damping)
        self.last_meas: VisionMeasurement | None = None
        self.filtered_center = np.zeros(3)
        self.filtered_velocity = np.zeros(3)
        self.baseline = np.zeros(3)
        self.last_update_time = 0.0

    def update_vision(self, meas: VisionMeasurement) -> None:
        """Accept one vision frame. A low-confidence frame is a dropout."""
        center = np.asarray(meas.center, dtype=float).reshape(3)
        if meas.confidence <= 0.4:
            return
        if self.last_meas is None:
            self.filtered_center = center.copy()
            self.filtered_velocity = np.zeros(3)
            self.baseline = center.copy()
        else:
            dt = meas.timestamp - self.last_meas.timestamp
            if dt <= 1e-5:
                return
            raw_velocity = (center - np.asarray(self.last_meas.center, dtype=float)) / dt
            self.filtered_velocity = 0.7 * self.filtered_velocity + 0.3 * raw_velocity
            self.filtered_center = 0.8 * self.filtered_center + 0.2 * center
            self.baseline = 0.99 * self.baseline + 0.01 * center
        self.last_meas = meas
        self.last_update_time = meas.timestamp

    def predict(self, current_time: float) -> tuple[Array, Array, float]:
        """Return centre (m), velocity (m/s), and an extra margin (m)."""
        if self.last_meas is None:
            return np.zeros(3), np.zeros(3), 0.0
        tau = float(np.clip(current_time - self.last_update_time, 0.0, 0.20))
        deviation, velocity = oscillator_step(
            self.filtered_center - self.baseline,
            self.filtered_velocity,
            tau,
            self.omega_b,
            self.damping,
        )
        # 0.5 mm at 100 ms of missing vision. The cap at 200 ms keeps this at 1 mm.
        uncertainty = 0.0005 * (tau / 0.10)
        return self.baseline + deviation, velocity, uncertainty


def oscillator_step(
    deviation: Array,
    velocity: Array,
    tau: float,
    omega: float,
    gamma: float,
) -> tuple[Array, Array]:
    """Exact step of ÿ + γ ẏ + ω² y = 0. Underdamped for a 0.3 Hz breath."""
    deviation = np.asarray(deviation, dtype=float).reshape(3)
    velocity = np.asarray(velocity, dtype=float).reshape(3)
    if tau <= 0.0:
        return deviation.copy(), velocity.copy()
    decay = 0.5 * gamma
    discriminant = omega**2 - decay**2
    if discriminant <= 1e-12:
        # Critically damped or overdamped is outside the breathing range used here.
        grown = np.exp(-gamma * tau)
        return deviation + velocity * tau * grown, velocity * grown
    rate = float(np.sqrt(discriminant))
    envelope = np.exp(-decay * tau)
    cosine = np.cos(rate * tau)
    sine = np.sin(rate * tau)
    amplitude = deviation
    slope = (velocity + decay * deviation) / rate
    next_deviation = envelope * (amplitude * cosine + slope * sine)
    next_velocity = -decay * next_deviation + envelope * (-amplitude * rate * sine + slope * rate * cosine)
    return next_deviation, next_velocity
