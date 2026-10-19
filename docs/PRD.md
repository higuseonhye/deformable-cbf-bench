# Product requirements

This document states what `deformable-cbf-bench` has to demonstrate. The repository is a synthetic kinematic benchmark. A service built on the same rules would still have to be validated on its own evidence. Measured values are the CPU run in [BENCHMARK_RESULTS.md](BENCHMARK_RESULTS.md).

## 1. Problem

A no-go sphere that stands still can be kept outside by a barrier that uses only the current centre. A sphere that moves into a stationary tool is a different failure. The half-space that drops tissue velocity stays quiet until the current barrier is already negative. One 100 Hz period of that delay spends the margin.

A second gap is the rate of the measurement. A vision sample may arrive at 30 Hz and then disappear for up to 100 ms. The command loop still asks for a centre, a velocity, and a shape every 10 ms.

A third gap is the surgical phase. A margin that is correct while approaching a vessel is the wrong margin at the moment the instrument is supposed to touch that vessel. The margin has to move. It must not jump inside one command.

## 2. Three mechanisms

### 2.1 Time-varying ellipsoid

For obstacle $`i`$, with centre $`p_i`$ and a symmetric shape matrix $`Q_i`$ that already contains the margin,

```math
h_i(x, t) = (x - p_i)^\top Q_i (x - p_i) - 1.
```

The safe set is $`h_i \ge 0`$. A sphere of anatomical radius $`r`$ and margin $`m`$ uses $`Q = I / (r + m)^2`$, so the zero level is the expanded sphere. Differentiating along $`\dot x = u`$ gives

```math
\dot h_i = 2(x - p_i)^\top Q_i (u - \dot p_i) + (x - p_i)^\top \dot Q_i (x - p_i).
```

Forward invariance asks for $`\dot h_i + \alpha h_i \ge 0`$. Rearranged for the command,

```math
A_i u \le b_i,
```

```math
A_i = -2 \bigl(Q_i (x - p_i)\bigr)^\top,
```

```math
b_i = \alpha h_i - 2(x - p_i)^\top Q_i \dot p_i + (x - p_i)^\top \dot Q_i (x - p_i).
```

The static call sets $`\dot p_i = 0`$ and $`\dot Q_i = 0`$. It still reacts once the current $`h_i`$ is negative. The nearest command is

```math
\begin{aligned}
\min_u \quad & \tfrac{1}{2}\|u - u_{\mathrm{nom}}\|^2 \\
\text{s.t.} \quad & A u \le b.
\end{aligned}
```

The scored plant is one sphere, a point tool, and a 10 ms Euler step. After the quadratic program, a repair adds speed along the outward normal, up to 40 mm/s, if the predicted barrier at the next step is still negative.

### 2.2 Vision gap

Between accepted frames the deviation $`y`$ from a slow baseline obeys

```math
\ddot y + \gamma \dot y + \omega_b^2 y = 0, \qquad \omega_b = 2\pi f_b.
```

The default breath is $`f_b = 0.3`$ Hz and $`\gamma = 0.95`$ s$`^{-1}`$. A frame with confidence at or below 0.4 does not reset the clock. The extra margin after a gap $`\tau`$ is

```math
\delta(\tau) = 0.5\,\mathrm{mm} \cdot \frac{\tau}{0.1\,\mathrm{s}}, \qquad 0 \le \tau \le 0.2\,\mathrm{s}.
```

This block is in metres. The hover and crossing episodes are in millimetres. The two are not added without a conversion.

### 2.3 Phase margin

Dissection uses a 1.0 mm margin. Hemostasis uses 0.3 mm. One control step blends

```math
d \leftarrow (1 - \lambda) d + \lambda d^\star, \qquad \lambda = \min(1, 5\Delta t).
```

At $`\Delta t = 10`$ ms, $`\lambda = 0.05`$. The force cap of the new phase is applied on that step. The position margin is not.

## 3. Users of this benchmark

- A batch path that reports hover clearance and crossing violations.
- A local board that replays the hover and one crossing, and draws the phase blend next to the dropout pad.

## 4. Functional requirements

| Item | Requirement |
| --- | --- |
| Plant | Point tool. One sphere, anatomical radius 5 mm, margin 1 mm |
| Hover | Nominal command 0. Sphere rises at 25 mm/s from the origin. Tool starts at z = 8 mm. 40 steps at 100 Hz |
| Static hover | Same integrator with tissue velocity removed from the half-space and from the one-step prediction |
| Crossing | Start (0, 18, 8) mm, target (0, −18, 8) mm. Centre z = 3 sin(2π · 0.3 t + φ) mm. Eight phases. 300 steps |
| Nominal command | Saturated attractor at 40 mm/s. Filtered runs add a tangent slide. `tv-replan` replaces the attractor for 0.5 s |
| Success | Clearance stays at least −0.001 mm and the final error is at most 4 mm |
| Solver | Active-set enumeration of $`A u \le b`$ for a 3-vector. Dual ascent is the fallback |
| Speed box | Separate metre-unit solver. $`\|u_j\| \le u_{\max}`$ stacked with the barrier. The box wins if the pair is empty |
| Predictor | Exact step of the linear oscillator. Dropout pad as in section 2.2 |
| Phase | Four named phases. Smooth margin and gain. Immediate force cap |

## 5. Out of scope

- A cleared device, a patient study, an endoscope recording, or a robot.
- A finite-element mesh in the loop. Vertex fitting is implemented and tested on a point cloud. The hover and crossing use a prescribed sphere.
- A trained flow field. The attractor is a saturated proportional command.
- Joint torque, cable backlash, and a vision model that names the phase.
- The Isaac Lab sketch under `deformable_cbf/envs`. Tests do not import it, and this run does not start a simulator.
- Treating the 0.3 Hz crossing as the episode that separates the two filters. On that bob the filters tie. The hover is the separating episode.

## 6. Acceptance

| Criterion | Target | Recorded run |
| --- | --- | --- |
| Static hover loses the margin | Clearance below 0 | −0.250 mm |
| Time-varying hover keeps the margin | Clearance at least 0 | +0.083 mm |
| Open-loop crossing loses the margin on some phases | At least one violation in 8 | 3 / 8, worst −0.994 mm |
| Filtered crossing keeps the margin | 0 violations | static, time-varying, and replan: 0 / 8 |
| Dropout pad at 100 ms | 0.5 mm, and a dim frame does not reset it | Unit test |
| Hemostasis margin does not jump in one 10 ms step | Strictly between 0.3 mm and 1.0 mm | 0.965 mm |
| Speed box | A 0.20 m/s command with no barrier is clipped to 0.05 m/s | Unit test |
