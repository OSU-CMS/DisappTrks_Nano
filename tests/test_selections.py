import pytest

ak = pytest.importorskip("awkward")

from disapptrks.selections import (
    build_lepton_veto_tag_probe_pairs,
    build_muon_veto_tag_probe_pairs,
)


def test_lepton_pairs_keep_probe_coordinates_for_fiducial_maps():
    tags = ak.Array(
        [
            [
                {
                    "pt": 50.0,
                    "eta": 0.0,
                    "phi": 0.0,
                    "charge": 1,
                }
            ]
        ]
    )
    probes = ak.Array(
        [
            [
                {
                    "pt": 45.0,
                    "eta": 1.25,
                    "phi": -2.5,
                    "charge": -1,
                    "dRMinElectron": 0.2,
                    "dRMinVetoElectron": 0.1,
                    "dRMinMuon": 0.3,
                    "dRMinTauHad": 0.4,
                    "dRMinJet": 0.6,
                    "caloEnergy": 5.0,
                    "missingOuterHits": 3,
                    "hp_trackerLayersWithMeasurement": 4,
                }
            ]
        ]
    )

    pairs = build_lepton_veto_tag_probe_pairs(
        tags,
        probes,
        tag_mass=0.000511,
        probe_mass=0.000511,
    )

    assert ak.to_list(pairs.probe_pt) == [[45.0]]
    assert ak.to_list(pairs.probe_eta) == [[1.25]]
    assert ak.to_list(pairs.probe_phi) == [[-2.5]]
    assert ak.to_list(pairs.probe_passElectronVeto) == [[True]]
    assert ak.to_list(pairs.probe_passVetoElectronVeto) == [[False]]


def test_muon_pairs_keep_loose_muon_veto_separate_from_generic_veto():
    tags = ak.Array(
        [
            [
                {
                    "pt": 50.0,
                    "eta": 0.0,
                    "phi": 0.0,
                    "charge": 1,
                }
            ]
        ]
    )
    probes = ak.Array(
        [
            [
                {
                    "pt": 45.0,
                    "eta": 1.25,
                    "phi": -2.5,
                    "charge": -1,
                    "dRMinMuon": 0.1,
                    "dRMinLooseMuon": 0.3,
                    "caloEnergy": 5.0,
                    "missingOuterHits": 3,
                    "hp_trackerLayersWithMeasurement": 4,
                }
            ]
        ]
    )

    pairs = build_muon_veto_tag_probe_pairs(tags, probes)

    assert ak.to_list(pairs.probe_passMuonVeto) == [[False]]
    assert ak.to_list(pairs.probe_passLooseMuonVeto) == [[True]]
