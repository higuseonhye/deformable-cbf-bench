# A kinematic check of a time-varying sphere barrier

This note rewrites one repository run. It was not peer reviewed. It does not claim a result on a surgical robot, on a finite-element organ, or on an endoscope recording.

## Abstract

On one CPU, a point tool was held still while a 5 mm sphere rose into it at 25 mm/s. A barrier that used the current centre and dropped the sphere velocity reached a clearance of −0.250 mm. The same quadratic program with the sphere velocity inside the half-space, plus a repair that refuses a negative barrier on the next 10 ms step, reached +0.083 mm. On a second episode the sphere bobbed by 3 mm at 0.3 Hz, peak speed 5.65 mm/s, and a saturated attractor drove the point through it. The attractor lost a 1 mm margin on 3 of 8 respiratory phases (worst clearance −0.994 mm). The static barrier and the time-varying barrier each lost it on 0 of 8 and reached the target on all 8. Replacing the attractor for 0.5 s did not change that count. Mean time around the filter was about 0.5 ms. A separate metre-unit layer integrates a damped breathing oscillator across a dropped vision frame, adding 0.5 mm of margin per 100 ms, and blends a phase change so a 1.0 mm margin becomes 0.965 mm on the first 10 ms step toward 0.3 mm rather than completing the jump.

## 1. Background

### 1.1 A barrier in front of a nominal command

A nominal velocity $`u_{\mathrm{nom}}`$ may be a hand-designed attractor or the output of a learned policy. The safety claim depends on the command that is integrated.

For a sphere that must not be entered, write the margin into the radius and set

```math
h(x, t) = \frac{\|x - p(t)\|^2}{(r + m)^2} - 1.
```

This is the spherical case of

```math
h_i(x, t) = (x - p_i)^\top Q_i (x - p_i) - 1,
```

with $`Q = I/(r+m)^2`$. The safe set $`\{h \ge 0\}`$ is the exterior. Along $`\dot x = u`$,

```math
\dot h = 2(x - p)^\top Q (u - \dot p) + (x - p)^\top \dot Q (x - p).
```

Invariance of the exterior asks for $`\dot h + \alpha h \ge 0`$. The quadratic program that stays nearest the nominal command is the single-integrator case of the control-barrier program in Ames, Xu, Grizzle, and Tabuada, IEEE Transactions on Automatic Control, 2017:

```math
\begin{aligned}
\min_u \quad & \tfrac{1}{2}\|u - u_{\mathrm{nom}}\|^2 \\
\text{s.t.} \quad & A u \le b.
\end{aligned}
```

If $`\dot p`$ is omitted, the half-space is the correct condition for a sphere that will still be where it is. It is the wrong condition for a sphere that is closing the gap during the control period.

### 1.2 A sample that is slower than the command

A vision measurement at 30 Hz leaves about three command ticks with no new centre. Smoke or an electrosurgical plume can extend that hole toward 100 ms. Holding the last position for the whole hole ignores the velocity already estimated from the previous pair of frames. Integrating that velocity without a bound lets a long dropout walk the obstacle into the tool. The bound used here is a linear oscillator at the breathing frequency, refused after 200 ms, with a margin that grows in proportion to the age of the last accepted frame.

### 1.3 A margin that depends on the phase

The same vessel is a region to stay clear of during dissection and a region the instrument must approach during hemostasis. A fixed 1 mm margin makes the second task infeasible. A switch that replaces 1.0 mm with 0.3 mm inside one 10 ms tick changes the admissible set discontinuously. The blend in section 3 moves 5% of the remaining gap on each tick.

## 2. Method

The scored plant is a point and one prescribed sphere. No mesh, no contact force, and no image is read. The command is a 3-vector. The quadratic program is solved by enumerating active sets. The step is then checked against the barrier value at $`t + 10`$ ms. A negative prediction adds outward speed up to 40 mm/s. The static arm is that procedure with the sphere velocity removed from the half-space and from the prediction. The open-loop arm is a saturated proportional attractor with neither the slide nor the filter.

The vision layer and the phase layer are not inside those episodes. They are checked by unit tests: a 100 ms dropout adds 0.5 mm and a low-confidence frame does not clear the timer; a zero-velocity offset of the oscillator moves back toward the baseline; the first hemostasis step changes the margin from 1.0 mm to 0.965 mm; a box constraint clips 0.20 m/s to 0.05 m/s.

## 3. One CPU run

Hover, 0.40 s, tissue speed 25 mm/s, nominal command 0:

| Filter | Minimum clearance |
| --- | ---: |
| Static | −0.250 mm |
| Time-varying | +0.083 mm |

−0.250 mm is one control period of the approach. The time-varying tool retreats and stays positive.

Crossing, eight phases, 3 mm bob at 0.3 Hz:

| Strategy | Margin lost | Success | Worst clearance |
| --- | ---: | ---: | ---: |
| Attractor | 3 / 8 | 62.5% | −0.994 mm |
| Static barrier | 0 / 8 | 100% | 0.848 mm |
| Time-varying barrier | 0 / 8 | 100% | 0.957 mm |
| Time-varying, with a 0.5 s replan | 0 / 8 | 100% | 0.957 mm |

Success means the margin holds and the final distance to the target is at most 4 mm. The peak tissue speed on this bob is 5.65 mm/s. The tangent slide is enough for both filters, which is why they tie. The hover is the episode that separates them.

## 4. Limits

The sphere velocity is prescribed. A tracker delay, a dropped vertex buffer, and a deformable finite-element step are not in the clearance numbers. The nominal policy is not a flow model. The replan direction is written by hand. Filter timing is one machine's mean, about half a millisecond, and is not a service-level percentile. The metre-unit predictor and the millimetre episodes do not share a state vector. Joint torque limits are not in the quadratic program; the box constraint is a bound on the Cartesian command.
