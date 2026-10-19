"""Local board for the moving-sphere barrier. Serves on 127.0.0.1:8506."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import streamlit as st

from deformable_cbf.bench import MARGIN_MM, RADIUS_MM
from deformable_cbf.traces import crossing_trace, dropout_margin_mm, hover_trace, phase_margin_mm

ROOT = Path(__file__).resolve().parent
REPORT_PATH = ROOT / "results" / "benchmark_report.json"

st.set_page_config(page_title="Moving-sphere barrier", layout="wide")


def _load_report() -> dict | None:
    if not REPORT_PATH.exists():
        return None
    raw = REPORT_PATH.read_text(encoding="utf-8").replace("NaN", "null")
    return json.loads(raw)


def _clearance_figure(static: dict, moving: dict):
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.axhline(0.0, color="#a8a29e", linewidth=1.0)
    ax.plot(static["t"], static["clearance"], color="#c2410c", linewidth=2.0, label="Static")
    ax.plot(moving["t"], moving["clearance"], color="#0f766e", linewidth=2.0, label="Time-varying")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Clearance (mm)")
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig


def _crossing_figure(attractor: dict, guided: dict):
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    radius = RADIUS_MM + MARGIN_MM
    for z, alpha in ((-3.0, 0.25), (0.0, 0.35), (3.0, 0.45)):
        circle = plt.Circle((0.0, z), radius, facecolor="#fecaca", edgecolor="#b91c1c", alpha=alpha, linewidth=0.8)
        ax.add_patch(circle)
    ax.plot(attractor["tool"][:, 1], attractor["tool"][:, 2], color="#64748b", linewidth=1.6, label="Attractor")
    ax.plot(guided["tool"][:, 1], guided["tool"][:, 2], color="#0f766e", linewidth=2.0, label="Time-varying")
    ax.scatter([18.0, -18.0], [8.0, 8.0], s=16, color="#0f172a", zorder=3)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("y (mm)")
    ax.set_ylabel("z (mm)")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    return fig


report = _load_report()
st.title("Moving-sphere barrier")
st.caption(
    "Synthetic point and a prescribed sphere. The cards are the committed CPU run. "
    "The sliders recompute that same integrator. Nothing here is a clinical device."
)

if report is None:
    st.warning("results/benchmark_report.json is missing. Run python -m deformable_cbf.bench.")
else:
    hover = {row["name"]: row for row in report["hover"]}
    crossing = {row["strategy"]: row for row in report["crossing"]}
    left, mid, right = st.columns(3)
    left.metric("Static hover clearance", f"{hover['hover-static']['min_clearance_mm']:.3f} mm")
    mid.metric("Time-varying hover clearance", f"{hover['hover-tv']['min_clearance_mm']:.3f} mm")
    right.metric("Attractor violations", f"{crossing['attractor']['violations']} / {crossing['attractor']['trials']}")

speed = st.slider("Hover tissue speed (mm/s)", min_value=5, max_value=40, value=25, step=1)
phase = st.slider("Crossing phase (rad)", min_value=0.0, max_value=6.28, value=0.0, step=0.01)

static = hover_trace(False, tissue_mm_s=float(speed))
moving = hover_trace(True, tissue_mm_s=float(speed))
attractor = crossing_trace("attractor", phase=phase)
guided = crossing_trace("tv", phase=phase)

column_a, column_b = st.columns(2)
with column_a:
    st.subheader("Hover")
    st.pyplot(_clearance_figure(static, moving))
    st.write(
        f"Live minimum clearance: static {float(static['clearance'].min()):.3f} mm, "
        f"time-varying {float(moving['clearance'].min()):.3f} mm. "
        f"The expanded radius is {RADIUS_MM + MARGIN_MM:.0f} mm."
    )
with column_b:
    st.subheader("Crossing")
    st.pyplot(_crossing_figure(attractor, guided))
    st.write(
        f"Live minimum clearance: attractor {float(attractor['clearance'].min()):.3f} mm, "
        f"time-varying {float(guided['clearance'].min()):.3f} mm."
    )

phase_t, phase_mm = phase_margin_mm()
gap_t, gap_mm = dropout_margin_mm()
schedule, schedule_ax = plt.subplots(figsize=(8.4, 2.8))
schedule_ax.plot(phase_t, phase_mm, color="#1d4ed8", linewidth=2.0, label="Margin after hemostasis")
schedule_ax.plot(gap_t, gap_mm, color="#7c3aed", linewidth=2.0, label="Pad after a dropped frame")
schedule_ax.set_xlabel("Time since the event (s)")
schedule_ax.set_ylabel("mm")
schedule_ax.legend(frameon=False)
schedule.tight_layout()
st.pyplot(schedule)
st.caption("The phase curve is in the metre package and drawn here in millimetres. The dropout pad is 0.5 mm per 100 ms, capped at 200 ms.")
