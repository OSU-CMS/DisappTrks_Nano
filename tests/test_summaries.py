from disapptrks.summaries import (
    cutflow_count,
    summarize_ss_subtracted_veto_probability,
    summarize_veto_probability,
)


def test_cutflow_count_sums_nested_nominal_counts():
    cutflow = {
        "muon_veto_zwindow": {
            "Run2024F": {
                "DATA_Muon": {"nominal": 6.0},
                "DATA_Other": {"nominal": 2.0},
            }
        }
    }

    assert cutflow_count(cutflow, "muon_veto_zwindow") == 8.0
    assert (
        cutflow_count(
            cutflow,
            "muon_veto_zwindow",
            dataset="Run2024F",
            sample="DATA_Muon",
        )
        == 6.0
    )


def test_summarize_veto_probability_uses_binomial_uncertainty():
    cutflow = {
        "den": {"dataset": {"sample": {"nominal": 100.0}}},
        "num": {"dataset": {"sample": {"nominal": 25.0}}},
    }

    summary = summarize_veto_probability(
        cutflow,
        denominator_name="den",
        numerator_name="num",
        dataset="dataset",
        sample="sample",
    )

    assert summary.denominator == 100.0
    assert summary.numerator == 25.0
    assert summary.probability == 0.25
    assert summary.uncertainty == (0.25 * 0.75 / 100.0) ** 0.5


def test_summarize_veto_probability_handles_empty_denominator():
    cutflow = {
        "den": {"dataset": {"sample": {"nominal": 0.0}}},
        "num": {"dataset": {"sample": {"nominal": 0.0}}},
    }

    summary = summarize_veto_probability(
        cutflow,
        denominator_name="den",
        numerator_name="num",
    )

    assert summary.probability == 0.0
    assert summary.uncertainty == 0.0


def test_summarize_ss_subtracted_veto_probability():
    cutflow = {
        "os_den": {"dataset": {"sample": {"nominal": 100.0}}},
        "os_num": {"dataset": {"sample": {"nominal": 30.0}}},
        "ss_den": {"dataset": {"sample": {"nominal": 20.0}}},
        "ss_num": {"dataset": {"sample": {"nominal": 10.0}}},
    }

    summary = summarize_ss_subtracted_veto_probability(
        cutflow,
        os_denominator_name="os_den",
        os_numerator_name="os_num",
        ss_denominator_name="ss_den",
        ss_numerator_name="ss_num",
        dataset="dataset",
        sample="sample",
    )

    assert summary.denominator == 80.0
    assert summary.numerator == 20.0
    assert summary.probability == 0.25
