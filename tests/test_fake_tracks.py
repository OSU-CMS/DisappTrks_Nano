import json

from disapptrks.fake_tracks import (
    estimate_fake_track_background,
    write_fake_track_latex,
)


def test_fake_track_estimate_from_cutflow_counts():
    cutflow = {
        "fake_basic3hits_d0_signal": {"dataset": {"sample": {"nominal": 20.0}}},
        "fake_basic3hits_d0_sideband": {"dataset": {"sample": {"nominal": 100.0}}},
        "fake_control_NLayers4": {"dataset": {"sample": {"nominal": 50.0}}},
    }

    estimate = estimate_fake_track_background(
        cutflow,
        layer="NLayers4",
        dataset="dataset",
        sample="sample",
    )

    assert estimate.transfer_factor.value == 0.2
    assert estimate.control.value == 50.0
    assert estimate.estimate.value == 10.0
    assert estimate.estimate.error > 0.0


def test_fake_track_estimate_applies_z_to_ll_normalization():
    counts = {
        "fake_basic3hits_d0_signal": 20.0,
        "fake_basic3hits_d0_sideband": 100.0,
        "fake_control_NLayers5": 50.0,
        "basic_yield": 2000.0,
        "ztoll_yield": 1000.0,
    }

    estimate = estimate_fake_track_background(
        counts,
        layer="NLayers5",
        basic_yield_category="basic_yield",
        z_to_ll_yield_category="ztoll_yield",
        source_is_cutflow=False,
    )

    assert estimate.normalization.value == 2.0
    assert estimate.control.value == 100.0
    assert estimate.estimate.value == 20.0
    assert estimate.p_fake_raw.value == 0.05


def test_write_fake_track_latex(tmp_path):
    counts = {
        "fake_basic3hits_d0_signal": 20.0,
        "fake_basic3hits_d0_sideband": 100.0,
        "fake_control_NLayers6plus": 50.0,
    }
    estimate = estimate_fake_track_background(
        counts,
        layer="NLayers6plus",
        source_is_cutflow=False,
    )

    path = tmp_path / "fake.tex"
    write_fake_track_latex([estimate], path, run_period="2023C")
    text = path.read_text()

    assert "run period" in text
    assert "2023C" in text
    assert "NLayers6plus" in text


def test_fake_track_estimate_json_serializable():
    counts = {
        "fake_basic3hits_d0_signal": 20.0,
        "fake_basic3hits_d0_sideband": 100.0,
        "fake_control_combinedBins": 50.0,
    }
    estimate = estimate_fake_track_background(
        counts,
        layer="combinedBins",
        source_is_cutflow=False,
    )

    json.dumps(estimate.as_dict())
