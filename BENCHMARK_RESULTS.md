# Benchmark results

One CPU run of `python -m deformable_cbf.bench`. The plant is a point mass and a prescribed sphere. Isaac Sim is not started.

## Hover, tissue rising at 25 mm/s, nominal command 0

| Filter | Min clearance (mm) | Margin lost |
| --- | ---: | --- |
| Static (ignore tissue velocity) | -0.250 | True |
| Time-varying | 0.083 | False |

## Crossing a 5 mm sphere, 3 mm / 0.3 Hz bob, 0.5 s occlusion

8 phases. Success means the margin holds and the final error is at most 4 mm.

| Strategy | Violations | Success | Worst clearance (mm) | Median final error (mm) | Mean filter time (us) |
| --- | ---: | ---: | ---: | ---: | ---: |
| attractor | 3/8 | 62.5% | -0.994 | 0.00 | 0.0 |
| static | 0/8 | 100.0% | 0.848 | 0.00 | 482.5 |
| tv | 0/8 | 100.0% | 0.957 | 0.00 | 528.2 |
| tv-replan | 0/8 | 100.0% | 0.957 | 0.01 | 541.2 |

Peak speed of the bob is 5.65 mm/s. On that crossing, the tangent slide is enough for both filters, so time-varying and static tie. The hover row is the separating case: the tissue approaches at 25 mm/s and the static filter reacts only after the current barrier is already negative.
