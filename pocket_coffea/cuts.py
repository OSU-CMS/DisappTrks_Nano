"""PocketCoffea event cuts for initial search-selection validation."""

from __future__ import annotations

import awkward as ak

from pocket_coffea.lib.cut_definition import Cut


def _has_disappearing_track(events, params, **kwargs):
    return events.nIsoTrackSearch >= params["minimum"]


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
