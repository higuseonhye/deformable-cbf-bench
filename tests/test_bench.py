from deformable_cbf.bench import crossing, run_suite


def test_open_loop_crossing_loses_the_margin_and_the_filter_does_not():
    lost = crossing("attractor", phase=0.0)
    kept = crossing("tv", phase=0.0)
    assert lost.violated
    assert not kept.violated


def test_suite_hover_and_tv_crossing():
    payload = run_suite(trials=4)
    hover = {row["name"]: row for row in payload["hover"]}
    assert hover["hover-static"]["violated"]
    assert not hover["hover-tv"]["violated"]
    by_name = {row["strategy"]: row for row in payload["crossing"]}
    assert by_name["attractor"]["violations"] > 0
    assert by_name["tv"]["violations"] == 0
    assert by_name["tv-replan"]["violations"] == 0
