import numpy as np

from disapptrks.fiducial import summarize_fiducial_map


def test_fiducial_hot_spot_is_identified():
    before = np.full((2, 2), 100.0)
    after = np.array([[1.0, 1.0], [1.0, 50.0]])
    edges = np.array([-1.0, 0.0, 1.0])

    summary = summarize_fiducial_map(
        before, after, edges, edges, threshold=1.0
    )

    assert len(summary.hot_spots) == 1
    assert summary.hot_spots[0].eta == 0.5
    assert summary.hot_spots[0].phi == 0.5
