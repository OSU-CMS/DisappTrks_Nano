from pathlib import Path

import numpy as np

from disapptrks.fake_tracks import Count
from disapptrks.lepton_backgrounds import (
    estimate_lepton_background,
    legacy_met_probabilities_from_outputs,
    write_lepton_background_json,
    write_lepton_background_latex,
)


class FakeAxis:
    def __init__(self, name, edges):
        self.name = name
        self.edges = np.asarray(edges, dtype=float)


class FakeHist:
    def __init__(self, counts, axes):
        self._counts = np.asarray(counts, dtype=float)
        self.axes = axes

    def values(self, flow=False):
        return self._counts


def test_estimate_lepton_background_uses_an_product():
    counts = {
        "electron_background_control_NLayers4": 100.0,
        "electron_background_offline_NLayers4": 25.0,
        "electron_background_trigger_NLayers4": 20.0,
    }
    pair_counts = {
        "NLayers4": {
            "den_os": 50.0,
            "num_os": 5.0,
            "den_ss": 10.0,
            "num_ss": 1.0,
        }
    }

    estimates = estimate_lepton_background(
        flavor=r"$e$",
        layers=["NLayers4"],
        pair_counts=pair_counts,
        counts=counts,
        control_category="electron_background_control_{layer}",
        poffline_numerator_category="electron_background_offline_{layer}",
        poffline_denominator_category="electron_background_control_{layer}",
        pmiss_numerator_category="electron_background_trigger_{layer}",
        pmiss_denominator_category="electron_background_offline_{layer}",
    )

    estimate = estimates[0]
    assert estimate.p_veto.value == 4.0 / 76.0
    assert estimate.p_offline.value == 0.25
    assert estimate.p_miss.value == 0.8
    assert estimate.p_trigger.value == 0.8
    assert estimate.estimate.value == 100.0 * (4.0 / 76.0) * 0.25 * 0.8


def test_estimate_lepton_background_applies_control_prescale_and_trigger_efficiency():
    estimates = estimate_lepton_background(
        flavor=r"$e$",
        layers=["NLayers4"],
        pair_counts={
            "NLayers4": {
                "den_os": 50.0,
                "num_os": 5.0,
                "den_ss": 10.0,
                "num_ss": 1.0,
            }
        },
        counts={
            "control_NLayers4": 100.0,
            "offline_NLayers4": 25.0,
            "trigger_NLayers4": 20.0,
        },
        control_category="control_{layer}",
        poffline_numerator_category="offline_{layer}",
        poffline_denominator_category="control_{layer}",
        pmiss_numerator_category="trigger_{layer}",
        pmiss_denominator_category="offline_{layer}",
        control_prescale=2.0,
        trigger_efficiency=Count(0.5, 0.0),
    )

    estimate = estimates[0]
    assert estimate.control_raw.value == 100.0
    assert estimate.control.value == 200.0
    assert estimate.trigger_efficiency.value == 0.5
    assert estimate.estimate.value == 200.0 * (4.0 / 76.0) * 0.25 * 0.8 / 0.5


def test_legacy_met_probabilities_integrate_trigger_turn_on():
    met_edges = [0.0, 120.0, 240.0]
    phi_edges = [0.0, 0.5, 1.0]
    output = {
        "variables": {
            "nElectronBackgroundMetMinusOnePt_NLayers4": {
                "sample": {
                    "dataset": FakeHist(
                        [10.0, 20.0],
                        [FakeAxis("met", met_edges)],
                    )
                }
            },
            "nElectronBackgroundMetMinusOnePtTrig_NLayers4": {
                "sample": {
                    "dataset": FakeHist(
                        [0.0, 10.0],
                        [FakeAxis("met", met_edges)],
                    )
                }
            },
            "nElectronBackgroundDeltaPhiMetJetLeadingVsMetMinusOnePt_NLayers4": {
                "sample": {
                    "dataset": FakeHist(
                        [[5.0, 5.0], [4.0, 6.0]],
                        [FakeAxis("met", met_edges), FakeAxis("dphi", phi_edges)],
                    )
                }
            },
        }
    }

    probabilities = legacy_met_probabilities_from_outputs(
        [output],
        prefix="Electron",
        layers=["NLayers4"],
        control_counts={"NLayers4": Count(40.0, 40.0)},
        dataset="dataset",
        sample="sample",
        met_cut=120.0,
        phi_cut=0.5,
    )

    poffline, pmiss = probabilities["NLayers4"]
    assert poffline.value == 6.0 / 40.0
    assert pmiss.value == 0.5


def test_write_lepton_background_outputs(tmp_path: Path):
    estimates = estimate_lepton_background(
        flavor=r"$\mu$",
        layers=["NLayers4"],
        pair_counts={
            "NLayers4": {
                "den_os": 10.0,
                "num_os": 1.0,
                "den_ss": 0.0,
                "num_ss": 0.0,
            }
        },
        counts={
            "control_NLayers4": 10.0,
            "offline_NLayers4": 5.0,
            "trigger_NLayers4": 4.0,
        },
        control_category="control_{layer}",
        poffline_numerator_category="offline_{layer}",
        poffline_denominator_category="control_{layer}",
        pmiss_numerator_category="trigger_{layer}",
        pmiss_denominator_category="offline_{layer}",
    )

    json_path = tmp_path / "lepton.json"
    tex_path = tmp_path / "lepton.tex"
    write_lepton_background_json(estimates, json_path)
    write_lepton_background_latex(estimates, tex_path, run_period="2023D")

    assert '"p_miss"' in json_path.read_text()
    text = tex_path.read_text()
    assert "\\epsilon_{\\mathrm{trig}}^{\\ell}" in text
    assert "P_{\\mathrm{offline}}" in text
    assert "P_{\\mathrm{miss}}" in text
    assert "2023D" in text
