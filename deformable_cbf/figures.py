"""Board figure for the recorded hover and one crossing phase.

Run `python -m deformable_cbf.figures` to rewrite docs/figures/dashboard.png.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from .bench import MARGIN_MM, RADIUS_MM
from .traces import crossing_trace, dropout_margin_mm, hover_trace, phase_margin_mm

ROOT = Path(__file__).resolve().parents[1]
FIGURE_PATH = ROOT / "docs" / "figures" / "dashboard.png"

STATIC = "#c2410c"
TV = "#0f766e"
OPEN = "#64748b"
SURFACE = "#b91c1c"
ZERO = "#a8a29e"


def render(path: Path = FIGURE_PATH) -> Path:
    static = hover_trace(False)
    moving = hover_trace(True)
    attractor = crossing_trace("attractor", phase=0.0)
    guided = crossing_trace("tv", phase=0.0)
    phase_t, phase_mm = phase_margin_mm()
    gap_t, gap_mm = dropout_margin_mm()

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 6.6))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.08, wspace=0.28, hspace=0.38)

    _hover_clearance(axes[0, 0], static, moving)
    _hover_height(axes[0, 1], static, moving)
    _crossing(axes[1, 0], attractor, guided)
    _schedules(axes[1, 1], phase_t, phase_mm, gap_t, gap_mm)

    fig.suptitle("Moving-sphere barrier on one CPU kinematic run", fontsize=13)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def _hover_clearance(ax, static: dict, moving: dict) -> None:
    ax.axhline(0.0, color=ZERO, linewidth=1.0)
    ax.plot(static["t"], static["clearance"], color=STATIC, linewidth=2.0, label="Static")
    ax.plot(moving["t"], moving["clearance"], color=TV, linewidth=2.0, label="Time-varying")
    ax.set_title("Hover clearance, tissue at 25 mm/s")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Clearance (mm)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(0.0, float(static["t"][-1]))


def _hover_height(ax, static: dict, moving: dict) -> None:
    ax.plot(static["t"], static["surface_z"], color=SURFACE, linewidth=1.6, label="Expanded surface")
    ax.plot(static["t"], static["tool_z"], color=STATIC, linewidth=2.0, label="Static tool")
    ax.plot(moving["t"], moving["tool_z"], color=TV, linewidth=2.0, label="Time-varying tool")
    ax.set_title("Same hover, height of the tool")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("z (mm)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(0.0, float(static["t"][-1]))


def _crossing(ax, attractor: dict, guided: dict) -> None:
    radius = RADIUS_MM + MARGIN_MM
    for z in (-3.0, 0.0, 3.0):
        ax.add_patch(Circle((0.0, z), radius, facecolor="#fecaca", edgecolor=SURFACE, alpha=0.35, linewidth=0.8))
    ax.plot(attractor["tool"][:, 1], attractor["tool"][:, 2], color=OPEN, linewidth=1.6, label="Attractor")
    ax.plot(guided["tool"][:, 1], guided["tool"][:, 2], color=TV, linewidth=2.0, label="Time-varying")
    ax.scatter([18.0, -18.0], [8.0, 8.0], s=18, color="#0f172a", zorder=3)
    ax.set_title("Crossing, phase 0, expanded sphere")
    ax.set_xlabel("y (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(frameon=False, fontsize=8, loc="lower left")


def _schedules(ax, phase_t, phase_mm, gap_t, gap_mm) -> None:
    ax.plot(phase_t, phase_mm, color="#1d4ed8", linewidth=2.0, label="Phase margin")
    ax.plot(gap_t, gap_mm, color="#7c3aed", linewidth=2.0, label="Dropout pad")
    ax.axhline(0.3, color="#1d4ed8", linewidth=0.7, linestyle="--")
    ax.set_title("Hemostasis blend and vision gap")
    ax.set_xlabel("Time since the event (s)")
    ax.set_ylabel("Millimetres")
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(0.0, 1.15)


def main() -> None:
    written = render()
    print(written)


if __name__ == "__main__":
    main()
