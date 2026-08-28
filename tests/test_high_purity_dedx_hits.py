import sys
from pathlib import Path

import awkward as ak

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pocket_coffea"))

from workflow import _dedx_hits_for_high_purity_tracks


def test_dedx_hits_are_linked_to_selected_isotrack_rows():
    events = ak.Array(
        [
            {
                "IsoTrackDeDxHit": [
                    {
                        "isoTrackIdx": 0,
                        "subdet": 1,
                        "layer": 2,
                        "isPixel": 1,
                        "passesStripShapeSelection": 1,
                        "pixelSize": 4,
                        "pixelSizeX": 2,
                        "pixelSizeY": 3,
                    },
                    {
                        "isoTrackIdx": 1,
                        "subdet": 5,
                        "layer": 3,
                        "isPixel": 0,
                        "passesStripShapeSelection": 0,
                        "pixelSize": -1,
                        "pixelSizeX": -1,
                        "pixelSizeY": -1,
                    },
                ]
            }
        ]
    )
    selected_tracks = ak.Array([[{"sourceIsoTrackIdx": 1}]])

    selected = _dedx_hits_for_high_purity_tracks(events, selected_tracks)

    assert ak.to_list(selected.isoTrackIdx) == [[1]]
    assert ak.to_list(selected.detectorLayer) == [[53]]
    assert ak.to_list(selected.pixelSize) == [[None]]
    assert ak.to_list(selected.stripPassesShapeSelection) == [[0]]
