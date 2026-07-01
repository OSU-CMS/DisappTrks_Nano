from pathlib import Path

from disapptrks.tables import (
    CountWithVariance,
    pveto_with_asymmetric_uncertainty,
    write_muon_cutflow_latex,
    write_muon_pveto_latex,
)


def test_pveto_with_asymmetric_uncertainty_uses_ss_subtraction():
    summary = pveto_with_asymmetric_uncertainty(
        den_os=CountWithVariance(100.0, 100.0),
        num_os=CountWithVariance(10.0, 10.0),
        den_ss=CountWithVariance(20.0, 20.0),
        num_ss=CountWithVariance(2.0, 2.0),
    )

    assert summary.denominator == 80.0
    assert summary.numerator == 8.0
    assert summary.central == 0.1
    assert summary.err_down > 0.0
    assert summary.err_up > 0.0


def test_pveto_with_negative_subtracted_numerator_quotes_upward_only():
    summary = pveto_with_asymmetric_uncertainty(
        den_os=CountWithVariance(100.0, 100.0),
        num_os=CountWithVariance(1.0, 1.0),
        den_ss=CountWithVariance(20.0, 20.0),
        num_ss=CountWithVariance(2.0, 2.0),
    )

    assert summary.numerator == -1.0
    assert summary.central == 0.0
    assert summary.err_down == 0.0
    assert summary.err_up > 0.0


def test_pveto_with_zero_subtracted_numerator_quotes_poisson_upper():
    summary = pveto_with_asymmetric_uncertainty(
        den_os=CountWithVariance(100.0, 100.0),
        num_os=CountWithVariance(0.0, 0.0),
        den_ss=CountWithVariance(0.0, 0.0),
        num_ss=CountWithVariance(0.0, 0.0),
    )

    assert summary.numerator == 0.0
    assert summary.denominator == 100.0
    assert summary.central == 0.0
    assert summary.err_down == 0.0
    assert summary.err_up > 0.0


def test_write_muon_latex_tables(tmp_path: Path):
    cutflow = {
        "inclusive": {"dataset": {"sample": {"nominal": 100.0}}},
        "muon_veto_tag": {"dataset": {"sample": {"nominal": 50.0}}},
        "muon_veto_probe": {"dataset": {"sample": {"nominal": 20.0}}},
        "muon_veto_pair": {"dataset": {"sample": {"nominal": 10.0}}},
        "muon_veto_pair_os_mass10": {"dataset": {"sample": {"nominal": 8.0}}},
        "muon_veto_zwindow": {"dataset": {"sample": {"nominal": 6.0}}},
        "muon_pveto_zwindow_pass": {"dataset": {"sample": {"nominal": 1.0}}},
        "muon_veto_ss_zwindow": {"dataset": {"sample": {"nominal": 2.0}}},
        "muon_pveto_ss_zwindow_pass": {"dataset": {"sample": {"nominal": 0.0}}},
        "muon_veto_zwindow_NLayers4": {"dataset": {"sample": {"nominal": 3.0}}},
        "muon_pveto_zwindow_pass_NLayers4": {"dataset": {"sample": {"nominal": 1.0}}},
        "muon_veto_ss_zwindow_NLayers4": {"dataset": {"sample": {"nominal": 1.0}}},
        "muon_pveto_ss_zwindow_pass_NLayers4": {"dataset": {"sample": {"nominal": 0.0}}},
    }

    cutflow_path = tmp_path / "cutflow.tex"
    pveto_path = tmp_path / "pveto.tex"
    write_muon_cutflow_latex(
        cutflow,
        cutflow_path,
        dataset="dataset",
        sample="sample",
    )
    summaries = write_muon_pveto_latex(
        cutflow,
        pveto_path,
        run_period="2024F",
        layers=["NLayers4", "combinedBins"],
        dataset="dataset",
        sample="sample",
    )

    assert "Cut & Events" in cutflow_path.read_text()
    assert "run period & flavor" in pveto_path.read_text()
    assert "2024F" in pveto_path.read_text()
    assert "N_{\\mathrm{layers}}=4" in pveto_path.read_text()
    assert summaries["combinedBins"].denominator == 4.0
    assert summaries["combinedBins"].numerator == 1.0
    assert summaries["NLayers4"].denominator == 2.0
    assert summaries["NLayers4"].numerator == 1.0
