from deformable_cbf.traces import hover_trace


def test_hover_trace_matches_the_recorded_clearance():
    static = hover_trace(False)
    moving = hover_trace(True)
    assert float(static["clearance"].min()) == -0.25
    assert abs(float(moving["clearance"].min()) - 0.08348762375258545) < 1e-9
