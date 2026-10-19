# Technical report

This report describes the CPU run scored in [BENCHMARK_RESULTS.md](BENCHMARK_RESULTS.md). It separates the barrier statement from the discrete step the runner integrates.

## 1. Plant

`bench.py` integrates a point at 100 Hz. The hover sphere has radius 5 mm and a 1 mm margin. It starts at the origin and moves at 25 mm/s along +z. The tool starts at z = 8 mm with a zero nominal command. Clearance is measured after the tool step and the sphere step:

```math
c = \|x - p\| - r - m.
```

The crossing sphere uses the same radius and margin. Its centre is

```math
p_z(t) = 3\sin(2\pi \cdot 0.3\, t + \varphi)\ \mathrm{mm},
```

so the peak speed is $`3 \cdot 2\pi \cdot 0.3 = 5.65`$ mm/s. The tool starts at y = 18 mm, z = 8 mm, and the target is y = −18 mm at the same height. Eight phases are equally spaced on $`[0, 2\pi)`$. An occlusion flag is true on steps 80 through 129. Only the `tv-replan` strategy changes its nominal command on that flag.

The nominal attractor is $`u = \mathrm{sat}(4(x^\star - x), 40\ \mathrm{mm/s})`$. Filtered strategies pass that command through a tangent slide when the squared barrier of the nearest obstacle is below 1 and the command points inward. The slide is part of the nominal command. It is not the certificate.

## 2. Barrier

### 2.1 Shape

`sphere_state` builds

```math
Q = \frac{1}{(r + m)^2} I, \qquad \dot Q = \frac{-2\dot r}{(r + m)^3} I.
```

The hover and the crossing have $`\dot r = 0`$, so $`\dot Q = 0`$. A general ellipsoid uses $`Q = R \operatorname{diag}(1/r_k^2) R^\top`$ after the radii have been expanded to contain the vertices and the margin has been added. That fit is in `tracker.py`. The scored episodes do not call it.

### 2.2 Half-space

With $`d = x - p`$,

```math
h = d^\top Q d - 1, \qquad A = -2(Q d)^\top,
```

```math
b = \alpha h - 2 d^\top Q \dot p + d^\top \dot Q d.
```

$`\alpha = 8`$ in the episodes. The static call zeros $`\dot p`$ and $`\dot Q`$ in both this half-space and the one-step prediction. The current centre is still the true centre, so a static filter that has already been entered will push back out. It does not push before the entry.

`solve_qp` enumerates active sets of $`A u \le b`$ because the command is a 3-vector and the row count on these episodes is 1. The stationarity condition used to accept a set is

```math
u - u_{\mathrm{nom}} + A^\top \lambda = 0, \qquad \lambda \ge 0.
```

If no set is feasible, or the row count exceeds 12, projected dual ascent runs on

```math
g(\lambda) = -\tfrac{1}{2}\lambda^\top A A^\top \lambda + \lambda^\top (A u_{\mathrm{nom}} - b)
```

with step $`1/\|AA^\top\|`$.

### 2.3 Discrete repair

The half-space is an instantaneous condition. The integrator then takes a finite step while the sphere also moves. `TVCBFFilter` predicts $`h`$ at $`x + \Delta t\, u`$ and at the Euler-advanced tissue state. While that value is negative it adds a step along $`Q d`$ of the most negative obstacle, of size at least $`0.2 u_{\max}`$, and caps the command at 40 mm/s. Up to 32 such additions are tried. Scaling a barely-safe outward command back toward zero is not used: when the tissue is the moving body, a shorter step is the less safe step.

On the hover, that repair is what keeps the time-varying tool on the safe side of a sphere closing at 25 mm/s. The static prediction does not include the 0.25 mm the sphere travels in one period, and the recorded clearance reaches −0.250 mm.

## 3. Crossing nominal

The open-loop strategy integrates the attractor with no slide and no filter. It lost the margin on 3 of 8 phases. The worst clearance was −0.994 mm. Every filtered strategy, including the static one, lost the margin on 0 of 8 and reached the target. The replan waypoint is a hand-written direction away from the nearest centre. On this bob it did not change the violation count. Mean filter time on the same run was 483 µs (static), 528 µs (time-varying), and 541 µs (replan). Those are wall-clock means on one machine, taken around the filter call, not a percentile from a server.

## 4. Metre-unit layer

`deformable_cbf/core` does not use the millimetre episodes.

The predictor stores an alpha-beta centre and velocity (weights 0.2 and 0.3 on the new sample) and a slower baseline (weight 0.01). The deviation from the baseline is propagated by the exact underdamped step of

```math
\ddot y + \gamma \dot y + \omega_b^2 y = 0.
```

With $`\gamma = 0.95`$ and $`\omega_b = 2\pi\cdot 0.3`$, the motion is underdamped. A zero-velocity offset therefore moves back toward the baseline. Extrapolation is refused past 200 ms. Confidence at or below 0.4 returns before the clock is touched, so the pad keeps growing. The pad is $`0.5`$ mm at 100 ms and $`1.0`$ mm at the cap.

The phase table is

| Phase | Margin | $`\alpha`$ | Force cap |
| --- | ---: | ---: | ---: |
| Approach | 1.5 mm | 20 | 1.0 N |
| Dissection | 1.0 mm | 25 | 2.5 N |
| Hemostasis | 0.3 mm | 40 | 4.0 N |
| Retraction | 0.8 mm | 15 | 3.0 N |

One 10 ms step from dissection toward hemostasis lands at

```math
0.95\cdot 1.0 + 0.05\cdot 0.3 = 0.965\ \mathrm{mm}.
```

The force cap on that step is already 4.0 N.

`BoxConstrainedTVCBFSolver` stacks $`|u_j| \le u_{\max}`$ with the barrier rows and calls the same active-set solver. If no point satisfies both, it clips the barrier-only solution to the box and marks the command saturated. A command of −0.20 m/s with an empty barrier is returned as −0.05 m/s.

## 5. What this run does not time

The vertex ellipsoid, the covariance pad, and the force row have unit tests and are not on the hover or crossing clock. The covariance pad adds $`\sqrt{\chi^2_3(0.999)}`$ times the largest standard deviation to every axis, with that square root equal to 4.033. It is a spectral buffer, not a tight chance constraint on $`h`$. The force row uses $`\dot F = K u`$ with $`K`$ in N/mm and the command in mm/s. No contact sensor produced the force in the test.

The Isaac Lab scene function is an import guarded by `ImportError`. It was not executed.
