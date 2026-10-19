# deformable-cbf-bench

Time-varying control barrier for ellipsoidal no-go zones that move.

The tool is a point with a velocity command, integrated at 100 Hz. Lengths are millimetres. For obstacle \(i\),

\[
h_i(x) = (x - p_i)^\top Q_i (x - p_i) - 1
\]

\(Q_i\) already includes the 1 mm margin. A sphere of anatomical radius \(r\) uses \(Q = I / (r + 1)^2\). The filter solves

\[
\min_u \tfrac12 \|u - u_{\text{nom}}\|^2 \quad \text{subject to} \quad \dot h_i + \alpha h_i \ge 0
\]

and then repairs the finite step. If the predicted \(h\) at \(t + \Delta t\) is still negative, speed is added along the outward normal, up to 40 mm/s. Slowing an outward command down is not a safe fallback when the tissue is what moves.

A static call is the same program with tissue velocity and \(\dot Q\) set to zero. It still reacts once the current \(h\) is negative. It does not act on the motion that will make \(h\) negative.

This repository does not start Isaac Sim, ORBIT, or a robot. The nominal command is a saturated attractor plus a tangent slide, not a trained flow field. Numbers in [BENCHMARK_RESULTS.md](BENCHMARK_RESULTS.md) are one CPU run of that kinematic check.

## Result that separates the two filters

The tool is told to hold still. A 5 mm sphere rises into it at 25 mm/s. Clearance is \(\|x-p\| - r - 1\).

| Filter | Min clearance | Margin lost |
| --- | ---: | --- |
| Static | −0.250 mm | yes |
| Time-varying | +0.083 mm | no |

−0.250 mm is one control period of that approach. The static half-space does not contain \(\dot p\), so the tool waits until the sphere has already crossed the margin.

On the slower crossing (3 mm bob at 0.3 Hz, peak tissue speed 5.65 mm/s) the two filters tie. Both keep the margin and reach the target on 8/8 phases. The unfiltered attractor loses the margin on 3/8. Swapping in the occlusion replan for 0.5 s does not change that outcome. The table is in [BENCHMARK_RESULTS.md](BENCHMARK_RESULTS.md).

Mean filter time on that run was about 0.5 ms, inside a 10 ms period. It is a wall-clock measurement of the active-set solve on this CPU, not a GPU kernel and not a P99 from a server.

## Layout

```
deformable_cbf/geometry.py   h, A, b and the one-step prediction
deformable_cbf/qp.py         active-set solve and a dual-ascent check
deformable_cbf/filter.py     projection plus outward repair
deformable_cbf/robust.py     pad each axis by sqrt(chi^2_3(0.999)) * sigma_max
deformable_cbf/force.py      one linear row, F_max^2 - ||F||^2, K in N/mm
deformable_cbf/tracker.py    containing ellipsoid of an (N, 3) vertex array
deformable_cbf/nominal.py    attractor, tangent slide, occlusion replan
deformable_cbf/bench.py      hover and crossing episodes
deformable_cbf/envs/scene_cfg.py
```

`scene_cfg.py` is an Isaac Lab 3 sketch (`volume_deformable_props`, `PhysxDeformableBodyMaterialCfg`). Nothing in the tests imports it. The tracker is the boundary the simulator would feed: nodal positions in, `EllipsoidState` out.

The covariance pad is a spectral buffer, not a tight chance constraint on \(h\). The force row uses \( \dot F = K u \) with the command expressed in the same frame as the force the tool applies to the tissue. It is not a measured organ law.

## Setup

Python 3.10+ and `requirements.txt`. NumPy is enough. A GPU is not used.

```bash
python -m pip install -r requirements.txt
pytest
python -m deformable_cbf.bench
```

## Not in this commit

Actuator torque limits, a margin that changes with surgical phase, and a predictor that bridges a 30 Hz vision frame to the 100 Hz command. Those are the next gaps. This filter assumes the ellipsoid state at the control tick is already known.
