import pytest

ak = pytest.importorskip("awkward")

from disapptrks.selections import (
    analysis_layer_mask,
    build_lepton_veto_tag_probe_pairs,
    build_muon_veto_tag_probe_pairs,
    fiducial_map_probe_track_mask,
    fake_track_no_d0_mask,
    hadronic_tau_control_object_mask,
    layer_mask,
    met_no_mu_minus_lepton,
    muon_veto_probe_track_mask,
    search_track_cutflow_masks,
    search_track_mask,
)


def test_analysis_layer_mask_requires_high_purity_only_for_four_layers():
    tracks = ak.Array(
        [[
            {"hp_trackerLayersWithMeasurement": 4, "isHighPurityTrack": False},
            {"hp_trackerLayersWithMeasurement": 4, "isHighPurityTrack": True},
            {"hp_trackerLayersWithMeasurement": 5, "isHighPurityTrack": False},
            {"hp_trackerLayersWithMeasurement": 6, "isHighPurityTrack": False},
        ]]
    )

    assert ak.to_list(analysis_layer_mask(tracks, "combinedBins")) == [
        [False, True, True, True]
    ]
    assert ak.to_list(analysis_layer_mask(tracks, "NLayers4")) == [
        [False, True, False, False]
    ]
    assert ak.to_list(analysis_layer_mask(tracks, "NLayers5")) == [
        [False, False, True, False]
    ]
    assert ak.to_list(analysis_layer_mask(tracks, "NLayers6plus")) == [
        [False, False, False, True]
    ]


def test_high_purity_study_can_retain_non_high_purity_four_layer_track():
    base = {
        "pt": 100.0, "eta": 0.8, "dxy": 0.1, "dz": 0.1,
        "inECALCrack": False, "inDTWheelGap": False,
        "inCSCTransition": False, "inTOBCrack": False,
        "isFiducialECALTrack": True, "hp_nValidPixelHits": 4,
        "hp_nValidHits": 4, "missingInnerHits": 0, "missingMiddleHits": 0,
        "missingOuterHits": 3, "pfRelIso03_chg": 0.01, "dRMinJet": 1.0,
        "dRMinElectron": 1.0, "dRMinMuon": 1.0, "dRMinTauHad": 1.0,
        "hp_trackerLayersWithMeasurement": 4, "caloEnergy": 1.0,
    }
    tracks = ak.Array([[
        {**base, "isHighPurityTrack": False},
        {**base, "isHighPurityTrack": True},
    ]])

    nominal = fake_track_no_d0_mask(tracks, layer="NLayers4")
    study = fake_track_no_d0_mask(
        tracks, layer="NLayers4", require_four_layer_high_purity=False
    )
    assert ak.to_list(nominal) == [[False, True]]
    assert ak.to_list(study) == [[True, True]]


def test_signal_track_comparison_changes_only_four_layer_high_purity():
    base = {
        "pt": 100.0, "eta": 0.8, "dxy": 0.01, "dz": 0.1,
        "inECALCrack": False, "inDTWheelGap": False,
        "inCSCTransition": False, "inTOBCrack": False,
        "isFiducialECALTrack": True, "hp_nValidPixelHits": 4,
        "hp_nValidHits": 4, "missingInnerHits": 0, "missingMiddleHits": 0,
        "missingOuterHits": 3, "pfRelIso03_chg": 0.01, "dRMinJet": 1.0,
        "dRMinElectron": 1.0, "dRMinMuon": 1.0, "dRMinTauHad": 1.0,
        "caloEnergy": 1.0,
    }
    tracks = ak.Array([[
        {
            **base,
            "hp_trackerLayersWithMeasurement": 4,
            "isHighPurityTrack": False,
        },
        {
            **base,
            "hp_trackerLayersWithMeasurement": 5,
            "isHighPurityTrack": False,
        },
    ]])

    nominal = search_track_mask(tracks)
    without_high_purity = search_track_mask(
        tracks,
        require_four_layer_high_purity=False,
    )

    assert ak.to_list(nominal) == [[False, True]]
    assert ak.to_list(without_high_purity) == [[True, True]]

    nominal_cutflow = search_track_cutflow_masks(tracks)
    comparison_cutflow = search_track_cutflow_masks(
        tracks,
        require_four_layer_high_purity=False,
    )
    assert ak.to_list(nominal_cutflow["track_layers4plus"]) == [[True, True]]
    assert ak.to_list(nominal_cutflow["track_highPurity4Layer"]) == [[False, True]]
    assert ak.to_list(comparison_cutflow["track_highPurity4Layer"]) == [[True, True]]


def test_signal_cartesian_axes_reproduce_post_high_purity_cutflows():
    base = {
        "pt": 100.0, "eta": 0.8, "dxy": 0.01, "dz": 0.1,
        "inECALCrack": False, "inDTWheelGap": False,
        "inCSCTransition": False, "inTOBCrack": False,
        "isFiducialECALTrack": True, "hp_nValidPixelHits": 4,
        "hp_nValidHits": 4, "missingInnerHits": 0, "missingMiddleHits": 0,
        "missingOuterHits": 3, "pfRelIso03_chg": 0.01, "dRMinJet": 1.0,
        "dRMinElectron": 1.0, "dRMinMuon": 1.0, "dRMinTauHad": 1.0,
        "caloEnergy": 1.0,
    }
    tracks = ak.Array([[
        {**base, "hp_trackerLayersWithMeasurement": 4, "isHighPurityTrack": False},
        {**base, "hp_trackerLayersWithMeasurement": 4, "isHighPurityTrack": True},
        {**base, "hp_trackerLayersWithMeasurement": 5, "isHighPurityTrack": False},
        {**base, "hp_trackerLayersWithMeasurement": 6, "isHighPurityTrack": False},
    ]])
    generic = search_track_cutflow_masks(
        tracks,
        layer="combinedBins",
        require_four_layer_high_purity=False,
    )
    post_high_purity_fields = list(generic)[
        list(generic).index("track_highPurity4Layer"):
    ]

    for layer in ("NLayers4", "NLayers5", "NLayers6plus", "combinedBins"):
        layer_axis = layer_mask(tracks, layer)
        for require_high_purity in (False, True):
            variant_axis = (
                (~layer_mask(tracks, "NLayers4") | tracks.isHighPurityTrack)
                if require_high_purity
                else ak.ones_like(tracks.pt, dtype=bool)
            )
            expected = search_track_cutflow_masks(
                tracks,
                layer=layer,
                require_four_layer_high_purity=require_high_purity,
            )
            for field in post_high_purity_fields:
                factored = generic[field] & layer_axis & variant_axis
                assert ak.to_list(factored) == ak.to_list(expected[field])


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


def test_met_no_mu_minus_tau_adds_visible_tau():
    events = ak.Array([{"MetNoMu": {"pt": 100.0, "phi": 0.0}}])
    taus = ak.Array([[{"pt": 50.0, "phi": 0.0}]])

    pt, phi = met_no_mu_minus_lepton(events, taus, flavor="tau")

    assert ak.to_list(pt) == pytest.approx([150.0])
    assert ak.to_list(phi) == pytest.approx([0.0])


def test_table27_tau_control_mask_uses_pt_eta_and_deeptau_raw_working_points():
    taus = ak.Array(
        [[
            {
                "pt": 60.0,
                "eta": 1.0,
                "idDecayModeNewDMs": True,
                "rawDeepTau2018v2p5VSjet": 0.90,
                "rawDeepTau2018v2p5VSe": 0.20,
                "rawDeepTau2018v2p5VSmu": 0.40,
            },
            {
                "pt": 49.0,
                "eta": 1.0,
                "idDecayModeNewDMs": True,
                "rawDeepTau2018v2p5VSjet": 0.90,
                "rawDeepTau2018v2p5VSe": 0.20,
                "rawDeepTau2018v2p5VSmu": 0.40,
            },
            {
                "pt": 60.0,
                "eta": 1.0,
                "idDecayModeNewDMs": True,
                "rawDeepTau2018v2p5VSjet": 0.80,
                "rawDeepTau2018v2p5VSe": 0.20,
                "rawDeepTau2018v2p5VSmu": 0.40,
            },
        ]]
    )

    assert ak.to_list(hadronic_tau_control_object_mask(taus)) == [
        [True, False, False]
    ]


def test_table27_tau_control_mask_accepts_osunano_wp_ordinals():
    taus = ak.Array(
        [[
            {
                "pt": 60.0,
                "eta": 1.0,
                "idDecayModeNewDMs": True,
                "idDeepTau2018v2p5VSjet": 6,
                "idDeepTau2018v2p5VSe": 1,
                "idDeepTau2018v2p5VSmu": 1,
            },
            {
                "pt": 60.0,
                "eta": 1.0,
                "idDecayModeNewDMs": True,
                "idDeepTau2018v2p5VSjet": 5,
                "idDeepTau2018v2p5VSe": 1,
                "idDeepTau2018v2p5VSmu": 1,
            },
        ]]
    )

    assert ak.to_list(hadronic_tau_control_object_mask(taus)) == [[True, False]]


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
