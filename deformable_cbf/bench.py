"""Kinematic episodes. No Isaac Sim, no learned policy, no recorded video.

Two checks:

1. Hover. The nominal command is zero while a sphere rises into the tool.
   A static barrier treats the sphere as fixed and stays put. The
   time-varying barrier retreats.

2. Crossing. A saturated attractor aims through a sphere that bobs on a
   0.3 Hz, 3 mm sinusoid. Occlusion lasts 0.5 s. REPLAN replaces the
   attractor during that window. The barrier still filters the result.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .filter import TVCBFFilter
from .geometry import sphere_clearance_mm, sphere_state
from .nominal import attractor, replan, slide

DT = 0.01
MARGIN_MM = 1.0
RADIUS_MM = 5.0
U_MAX = 40.0
ALPHA = 8.0


@dataclass(frozen=True)
class EpisodeResult:
    name: str
    min_clearance_mm: float
    violated: bool
    final_error_mm: float
    reached: bool
    interventions: int
    mean_filter_us: float


def hover(predict_motion: bool, steps: int = 40, tissue_mm_s: float = 25.0) -> EpisodeResult:
    """Tool wants to hold. Tissue velocity is +z, into the tool."""
    filt = TVCBFFilter(alpha=ALPHA, dt=DT, u_max=U_MAX)
    position = np.array([0.0, 0.0, 8.0])
    center = np.array([0.0, 0.0, 0.0])
    velocity = np.array([0.0, 0.0, tissue_mm_s])
    clearances = []
    interventions = 0
    for _ in range(steps):
        state = sphere_state(center.reshape(1, 3), velocity.reshape(1, 3), np.array([RADIUS_MM]), MARGIN_MM)
        result = filt.apply(position, np.zeros(3), state, predict_motion=predict_motion)
        if result.intervened:
            interventions += 1
        position = position + DT * result.velocity
        center = center + DT * velocity
        clearances.append(sphere_clearance_mm(position, center, RADIUS_MM, MARGIN_MM))
    minimum = float(np.min(clearances))
    return EpisodeResult(
        name="hover-tv" if predict_motion else "hover-static",
        min_clearance_mm=minimum,
        violated=bool(minimum < -1e-3),
        final_error_mm=float("nan"),
        reached=False,
        interventions=interventions,
        mean_filter_us=0.0,
    )


def crossing(strategy: str, phase: float, steps: int = 300, occlusion: tuple[int, int] = (80, 130)) -> EpisodeResult:
    filt = TVCBFFilter(alpha=ALPHA, dt=DT, u_max=U_MAX)
    position = np.array([0.0, 18.0, 8.0])
    target = np.array([0.0, -18.0, 8.0])
    clearances = []
    interventions = 0
    times = []
    for step in range(steps):
        t = step * DT
        z = 3.0 * np.sin(2.0 * np.pi * 0.3 * t + phase)
        vz = 3.0 * (2.0 * np.pi * 0.3) * np.cos(2.0 * np.pi * 0.3 * t + phase)
        center = np.array([0.0, 0.0, z])
        state = sphere_state(center.reshape(1, 3), np.array([[0.0, 0.0, vz]]), np.array([RADIUS_MM]), MARGIN_MM)
        occluded = occlusion[0] <= step < occlusion[1]
        command = _nominal(strategy, position, target, state, occluded)
        if strategy == "attractor":
            applied = command
        else:
            started = _ticks()
            result = filt.apply(position, command, state, predict_motion=strategy != "static")
            times.append(_ticks() - started)
            applied = result.velocity
            if result.intervened:
                interventions += 1
        position = position + DT * applied
        t_next = (step + 1) * DT
        z_next = 3.0 * np.sin(2.0 * np.pi * 0.3 * t_next + phase)
        clearances.append(sphere_clearance_mm(position, np.array([0.0, 0.0, z_next]), RADIUS_MM, MARGIN_MM))
    minimum = float(np.min(clearances))
    final_error = float(np.linalg.norm(position - target))
    return EpisodeResult(
        name=strategy,
        min_clearance_mm=minimum,
        violated=bool(minimum < -1e-3),
        final_error_mm=final_error,
        reached=bool(final_error <= 4.0 and minimum >= -1e-3),
        interventions=interventions,
        mean_filter_us=float(np.mean(times) * 1e6) if times else 0.0,
    )


def run_suite(trials: int = 8, seed: int = 0) -> dict:
    phases = np.linspace(0.0, 2.0 * np.pi, trials, endpoint=False)
    hover_rows = [hover(False), hover(True)]
    crossing_rows = []
    for strategy in ("attractor", "static", "tv", "tv-replan"):
        runs = [crossing(strategy, float(phase)) for phase in phases]
        crossing_rows.append(_aggregate(strategy, runs))
    return {
        "dt_s": DT,
        "margin_mm": MARGIN_MM,
        "radius_mm": RADIUS_MM,
        "u_max_mm_s": U_MAX,
        "alpha": ALPHA,
        "trials": trials,
        "seed": seed,
        "hover": [asdict(row) for row in hover_rows],
        "crossing": crossing_rows,
        "notes": [
            "Kinematic single integrator. The sphere centre is prescribed, not simulated by FEM.",
            "The nominal command is a saturated attractor, a tangent slide, and a hand-written replan.",
            "Clearance is ||x - p|| - radius - margin after the tissue's Euler step.",
        ],
    }


def write_report(payload: dict, root: Path) -> None:
    result_dir = root / "results"
    result_dir.mkdir(exist_ok=True)
    (result_dir / "benchmark_report.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (root / "BENCHMARK_RESULTS.md").write_text(_markdown(payload), encoding="utf-8")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    payload = run_suite()
    write_report(payload, root)
    print(_markdown(payload))


def _nominal(strategy: str, position: np.ndarray, target: np.ndarray, state, occluded: bool) -> np.ndarray:
    if strategy == "tv-replan" and occluded:
        return replan(position, state, U_MAX)
    command = attractor(position, target, gain=4.0, u_max=U_MAX)
    if strategy == "attractor":
        return command
    return slide(position, command, state, target, U_MAX, alpha=ALPHA)


def _aggregate(strategy: str, runs: list[EpisodeResult]) -> dict:
    count = len(runs)
    violations = sum(row.violated for row in runs)
    reached = sum(row.reached for row in runs)
    return {
        "strategy": strategy,
        "trials": count,
        "violations": violations,
        "violation_rate": violations / count,
        "reached": reached,
        "success_rate": reached / count,
        "min_clearance_mm": float(np.min([row.min_clearance_mm for row in runs])),
        "median_final_error_mm": float(np.median([row.final_error_mm for row in runs])),
        "mean_filter_us": float(np.mean([row.mean_filter_us for row in runs])),
    }


def _markdown(payload: dict) -> str:
    hover = {row["name"]: row for row in payload["hover"]}
    lines = [
        "# Benchmark results",
        "",
        "One CPU run of `python -m deformable_cbf.bench`. The plant is a point mass and a prescribed sphere. Isaac Sim is not started.",
        "",
        "## Hover, tissue rising at 25 mm/s, nominal command 0",
        "",
        "| Filter | Min clearance (mm) | Margin lost |",
        "| --- | ---: | --- |",
        f"| Static (ignore tissue velocity) | {hover['hover-static']['min_clearance_mm']:.3f} | {hover['hover-static']['violated']} |",
        f"| Time-varying | {hover['hover-tv']['min_clearance_mm']:.3f} | {hover['hover-tv']['violated']} |",
        "",
        "## Crossing a 5 mm sphere, 3 mm / 0.3 Hz bob, 0.5 s occlusion",
        "",
        f"{payload['trials']} phases. Success means the margin holds and the final error is at most 4 mm.",
        "",
        "| Strategy | Violations | Success | Worst clearance (mm) | Median final error (mm) | Mean filter time (us) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in payload["crossing"]:
        lines.append(
            f"| {row['strategy']} | {row['violations']}/{row['trials']} | {row['success_rate'] * 100:.1f}% | {row['min_clearance_mm']:.3f} | {row['median_final_error_mm']:.2f} | {row['mean_filter_us']:.1f} |"
        )
    lines.extend(
        [
            "",
            "Peak speed of the bob is 5.65 mm/s. On that crossing, the tangent slide is enough for both filters, so time-varying and static tie. The hover row is the separating case: the tissue approaches at 25 mm/s and the static filter reacts only after the current barrier is already negative.",
            "",
        ]
    )
    return "\n".join(lines)


def _ticks() -> float:
    return __import__("time").perf_counter()
