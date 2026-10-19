# Benchmark results

This is the CPU run also stored at the repository root and in `results/benchmark_report.json`. The plant is a point mass and a prescribed sphere. A simulator is not started. Filter times are wall-clock means around the filter call on that run.

## Hover, tissue rising at 25 mm/s, nominal command 0

| Filter | Min clearance (mm) | Margin lost |
| --- | ---: | --- |
| Static (tissue velocity omitted) | −0.250 | yes |
| Time-varying | +0.083 | no |

## Crossing a 5 mm sphere, 3 mm at 0.3 Hz, 0.5 s occlusion flag

Eight phases. Success means the margin holds and the final error is at most 4 mm.

| Strategy | Violations | Success | Worst clearance (mm) | Median final error (mm) | Mean filter time (µs) |
| --- | ---: | ---: | ---: | ---: | ---: |
| attractor | 3 / 8 | 62.5% | −0.994 | 0.00 | — |
| static | 0 / 8 | 100% | 0.848 | 0.00 | 483 |
| tv | 0 / 8 | 100% | 0.957 | 0.00 | 528 |
| tv-replan | 0 / 8 | 100% | 0.957 | 0.01 | 541 |

Peak speed of the bob is 5.65 mm/s. On that crossing the tangent slide is enough for both filters. The hover row is the separating case.
