from disapptrks.datasets import build_dataset_definition, root_files_from_lines


def test_root_files_from_lines_filters_comments_and_non_root_files():
    files = root_files_from_lines(
        [
            "",
            "# comment",
            "root://cmseos.fnal.gov//store/a/file2.root",
            "root://cmseos.fnal.gov//store/a/not_root.txt",
            "root://cmseos.fnal.gov//store/a/file1.root",
        ]
    )

    assert files == [
        "root://cmseos.fnal.gov//store/a/file1.root",
        "root://cmseos.fnal.gov//store/a/file2.root",
    ]


def test_build_dataset_definition_includes_required_metadata():
    dataset = build_dataset_definition(
        dataset_name="Run2024G_Muon_OSUNano_EOS",
        files=["root://cmseos.fnal.gov//store/a/file.root"],
        sample="DATA_Muon",
        year="2024",
        era="G",
        primary_dataset="Muon",
        nevents="0",
    )

    entry = dataset["Run2024G_Muon_OSUNano_EOS"]
    assert entry["files"] == ["root://cmseos.fnal.gov//store/a/file.root"]
    assert entry["metadata"]["sample"] == "DATA_Muon"
    assert entry["metadata"]["year"] == "2024"
    assert entry["metadata"]["era"] == "G"
    assert entry["metadata"]["nevents"] == "0"
    assert entry["metadata"]["isMC"] == "False"
