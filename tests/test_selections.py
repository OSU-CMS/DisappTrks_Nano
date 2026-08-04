import pytest

ak = pytest.importorskip("awkward")

from disapptrks.selections import (
    build_lepton_veto_tag_probe_pairs,
    build_muon_veto_tag_probe_pairs,
    fiducial_map_probe_track_mask,
    met_no_mu_minus_lepton,
    muon_veto_probe_track_mask,
)


def test_met_no_mu_minus_muon_does_not_add_muon_twice():
    events = ak.Array([{"MetNoMu": {"pt": 100.0, "phi": 0.0}}])
    muons = ak.Array([[{"pt": 40.0, "phi": 0.0}]])

    pt, phi = met_no_mu_minus_lepton(events, muons, flavor="muon")

    assert ak.to_list(pt) == pytest.approx([100.0])
    assert ak.to_list(phi) == pytest.approx([0.0])


def test_met_no_mu_minus_electron_adds_visible_tag():
    events = ak.Array([{"MetNoMu": {"pt": 100.0, "phi": 0.0}}])
    electrons = ak.Array([[{"pt": 40.0, "phi": 0.0}]])

    pt, phi = met_no_mu_minus_lepton(events, electrons, flavor="electron")

    assert ak.to_list(pt) == pytest.approx([140.0])
    assert ak.to_list(phi) == pytest.approx([0.0])


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


def test_fiducial_map_probe_uses_legacy_old_hit_cuts():
    tracks = ak.Array(
        [
            [
                {
                    "pt": 45.0,
                    "eta": 0.6,
                    "phi": 1.0,
                    "inECALCrack": False,
                    "inDTWheelGap": False,
                    "inCSCTransition": False,
                    "inTOBCrack": False,
                    "isFiducialECALTrack": True,
                    "hp_nValidPixelHits": 3,
                    "hp_nValidHits": 7,
                    "missingInnerHits": 0,
                    "missingMiddleHits": 0,
                    "pfRelIso03_chg": 0.01,
                    "dxy": 0.01,
                    "dz": 0.1,
                    "hp_trackerLayersWithMeasurement": 4,
                    "dRMinJet": 0.6,
                    "dRMinElectron": 0.2,
                    "dRMinMuon": 0.2,
                    "dRMinTauHad": 0.2,
                    "caloEnergy": 5.0,
                }
            ]
        ]
    )

    assert ak.to_list(fiducial_map_probe_track_mask(tracks, flavor="muon")) == [
        [True]
    ]
    assert ak.to_list(muon_veto_probe_track_mask(tracks)) == [[False]]
