from pathlib import Path

import pytest

from disapptrks.cli import _period_float_values, _standard_coffea_files


def test_standard_coffea_files_prefers_merged_output(tmp_path: Path):
    merged = tmp_path / "output_all.coffea"
    shard = tmp_path / "output_job_0.coffea"
    merged.touch()
    shard.touch()

    assert _standard_coffea_files(tmp_path) == [merged]


def test_standard_coffea_files_uses_sorted_job_shards(tmp_path: Path):
    second = tmp_path / "output_job_2.coffea"
    first = tmp_path / "output_job_1.coffea"
    second.touch()
    first.touch()

    assert _standard_coffea_files(tmp_path) == [first, second]


def test_standard_coffea_files_rejects_empty_directory(tmp_path: Path):
    with pytest.raises(SystemExit, match="no top-level .coffea files"):
        _standard_coffea_files(tmp_path)


def test_standard_coffea_files_rejects_ambiguous_outputs(tmp_path: Path):
    (tmp_path / "output_category_a.coffea").touch()
    (tmp_path / "output_category_b.coffea").touch()

    with pytest.raises(SystemExit, match="ambiguous .coffea outputs"):
        _standard_coffea_files(tmp_path)


def test_period_float_values_accepts_single_bare_value():
    assert _period_float_values(["0.9"], ["2022CD"], "--value") == {
        "2022CD": 0.9
    }


def test_period_float_values_maps_multiple_periods():
    assert _period_float_values(
        ["2022CD=0.9", "2022EFG=0.91"],
        ["2022CD", "2022EFG"],
        "--value",
    ) == {"2022CD": 0.9, "2022EFG": 0.91}


def test_period_float_values_rejects_missing_period():
    with pytest.raises(SystemExit, match="missing for: 2022EFG"):
        _period_float_values(
            ["2022CD=0.9"], ["2022CD", "2022EFG"], "--value"
        )
