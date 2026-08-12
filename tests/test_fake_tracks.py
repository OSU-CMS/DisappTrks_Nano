import json

from disapptrks.fake_tracks import (
    estimate_fake_track_background,
    fixed_an_transfer_factor_fit,
    write_combined_fake_track_table34_latex,
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


def test_fixed_an_transfer_factor_accepts_run_period_aliases():
    fit = fixed_an_transfer_factor_fit("2022 CD", "zmumu")

    assert fit.histogram == "fixed:2022CD:zmumu"
    assert fit.transfer_factor.value == 0.10
    assert fit.transfer_factor.error == 0.06


def test_combined_table34_contains_each_run_period(tmp_path):
    period_paths = {}
    for period in ("2022CD", "2022EFG"):
        json_paths = []
        for control_region in (r"$Z\to\mu\mu$", r"$Z\to ee$"):
            path = tmp_path / f"{period}_{len(json_paths)}.json"
            path.write_text(
                json.dumps(
                    {
                        "estimates": [
                            {
                                "layer": "NLayers4",
                                "control_region": control_region,
                                "fake_probability": {"value": 0.001, "variance": 1e-8},
                                "fake_yield": {"value": 2.0, "variance": 0.25},
                            }
                        ]
                    }
                )
            )
            json_paths.append(path)
        period_paths[period] = json_paths

    output = tmp_path / "combined.tex"
    write_combined_fake_track_table34_latex(period_paths, output)

    text = output.read_text()
    assert "2022CD" in text
    assert "2022EFG" in text
