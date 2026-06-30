"""PocketCoffea event cuts for initial search-selection validation."""

from __future__ import annotations

import awkward as ak

from pocket_coffea.lib.cut_definition import Cut


EVENT_DIAGNOSTIC_FIELDS = [
    "event_metNoMu120",
    "event_leadingJet110",
    "event_jetMetDphi0p5",
    "event_dijetDphi2p5",
]

TRACK_DIAGNOSTIC_FIELDS = [
    "track_pt55",
    "track_eta2p1",
    "track_noECALCrack",
    "track_noDTWheelGap",
    "track_noCSCTransition",
    "track_noTOBCrack",
    "track_fiducialECAL",
    "track_pixelHits4",
    "track_validHits4",
    "track_noMissingInner",
    "track_noMissingMiddle",
    "track_chargedIso0p05",
    "track_dxy0p02",
    "track_dz0p5",
    "track_dRJet0p5",
    "track_layers4plus",
    "track_calo10",
    "track_missingOuter3",
    "track_electronVeto",
    "track_muonVeto",
    "track_tauVeto",
]

COMBINED_DIAGNOSTIC_FIELDS = [
    f"eventKinematics_{field}" for field in TRACK_DIAGNOSTIC_FIELDS
]

SEARCH_DIAGNOSTIC_FIELDS = (
    EVENT_DIAGNOSTIC_FIELDS + TRACK_DIAGNOSTIC_FIELDS + COMBINED_DIAGNOSTIC_FIELDS
)


def _has_disappearing_track(events, params, **kwargs):
    return events.nIsoTrackSearch >= params["minimum"]


def _has_count(events, params, **kwargs):
    return events[params["field"]] >= params["minimum"]


def _search_diagnostic(events, params, **kwargs):
    return events.SearchDiag[params["field"]]


def _search_kinematics(events, params, **kwargs):
    event = events.AnalysisEvent
    return (
        (event.METNoMu_pt >= params["met_min"])
        & (event.leadingJet_pt > params["jet_pt_min"])
        & (event.leadingJetMETNoMuDeltaPhi >= params["jet_met_dphi_min"])
        & ((event.dijetMaxDeltaPhi < 0.0) | (event.dijetMaxDeltaPhi < params["dijet_dphi_max"]))
    )


has_disappearing_track = Cut(
    name="has_disappearing_track",
    params={"minimum": 1},
    function=_has_disappearing_track,
)

has_muon_tag = Cut(
    name="has_muon_tag",
    params={"field": "nMuonTag", "minimum": 1},
    function=_has_count,
)

has_muon_veto_probe_track = Cut(
    name="has_muon_veto_probe_track",
    params={"field": "nMuonVetoProbeTrack", "minimum": 1},
    function=_has_count,
)

has_muon_veto_tag_probe_pair = Cut(
    name="has_muon_veto_tag_probe_pair",
    params={"field": "nMuonVetoTagProbePair", "minimum": 1},
    function=_has_count,
)

has_muon_veto_zwindow_pair = Cut(
    name="has_muon_veto_zwindow_pair",
    params={"field": "nMuonVetoTagProbePairZWindow", "minimum": 1},
    function=_has_count,
)

has_muon_veto_zwindow_pass_pair = Cut(
    name="has_muon_veto_zwindow_pass_pair",
    params={"field": "nMuonVetoTagProbePairZWindowPass", "minimum": 1},
    function=_has_count,
)

search_kinematics = Cut(
    name="search_kinematics",
    params={
        "met_min": 120.0,
        "jet_pt_min": 110.0,
        "jet_met_dphi_min": 0.5,
        "dijet_dphi_max": 2.5,
    },
    function=_search_kinematics,
)

search_diagnostic_cuts = {
    field: Cut(
        name=field,
        params={"field": field},
        function=_search_diagnostic,
    )
    for field in SEARCH_DIAGNOSTIC_FIELDS
}
