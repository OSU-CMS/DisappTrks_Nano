"""PocketCoffea workflow for the custom disappearing-track NanoAOD."""

from __future__ import annotations

import awkward as ak

from pocket_coffea.workflows.base import BaseProcessorABC

from disapptrks.selections import (
    add_event_derived_fields,
    add_isotrack_derived_fields,
    base_probe_track_mask,
    build_muon_veto_tag_probe_pairs,
    mass10_muon_probe_pair_mask,
    muon_tag_mask,
    muon_veto_probe_track_mask,
    muon_veto_pair_fail_mask,
    muon_veto_pair_pass_mask,
    os_muon_probe_pair_mask,
    search_event_cutflow_masks,
    search_track_cutflow_masks,
    search_track_mask,
    z_window_muon_probe_pair_mask,
)


class DisappTrksProcessor(BaseProcessorABC):
    def apply_object_preselection(self, variation):
        self.events["IsoTrack"] = add_isotrack_derived_fields(self.events)
        self.events["AnalysisEvent"] = add_event_derived_fields(self.events)
        self.events["MuonTag"] = self.events.Muon[muon_tag_mask(self.events.Muon)]
        self.events["IsoTrackProbe"] = self.events.IsoTrack[
            base_probe_track_mask(self.events.IsoTrack)
        ]
        self.events["MuonVetoProbeTrack"] = self.events.IsoTrack[
            muon_veto_probe_track_mask(self.events.IsoTrack)
        ]
        muon_veto_pairs = build_muon_veto_tag_probe_pairs(
            self.events.MuonTag, self.events.MuonVetoProbeTrack
        )
        self.events["MuonVetoTagProbePair"] = muon_veto_pairs
        self.events["MuonVetoTagProbePairOS"] = muon_veto_pairs[
            os_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairOSMass10"] = muon_veto_pairs[
            mass10_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairZWindow"] = muon_veto_pairs[
            z_window_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairZWindowPass"] = muon_veto_pairs[
            z_window_muon_probe_pair_mask(muon_veto_pairs)
            & muon_veto_pair_pass_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairZWindowFail"] = muon_veto_pairs[
            z_window_muon_probe_pair_mask(muon_veto_pairs)
            & muon_veto_pair_fail_mask(muon_veto_pairs)
        ]
        search_diagnostic_masks = search_track_cutflow_masks(self.events.IsoTrack)
        self.events["IsoTrackSearchPreMissingOuter"] = self.events.IsoTrack[
            search_diagnostic_masks["track_calo10"]
        ]
        self.events["IsoTrackSearchPreLeptonVeto"] = self.events.IsoTrack[
            search_diagnostic_masks["track_missingOuter3"]
        ]
        self.events["IsoTrackSearch"] = self.events.IsoTrack[
            search_track_mask(self.events.IsoTrack)
        ]

    def count_objects(self, variation):
        self.events["nIsoTrack"] = ak.num(self.events.IsoTrack)
        self.events["nMuonTag"] = ak.num(self.events.MuonTag)
        self.events["nIsoTrackProbe"] = ak.num(self.events.IsoTrackProbe)
        self.events["nMuonVetoProbeTrack"] = ak.num(self.events.MuonVetoProbeTrack)
        self.events["nMuonVetoTagProbePair"] = ak.num(
            self.events.MuonVetoTagProbePair
        )
        self.events["nMuonVetoTagProbePairOS"] = ak.num(
            self.events.MuonVetoTagProbePairOS
        )
        self.events["nMuonVetoTagProbePairOSMass10"] = ak.num(
            self.events.MuonVetoTagProbePairOSMass10
        )
        self.events["nMuonVetoTagProbePairZWindow"] = ak.num(
            self.events.MuonVetoTagProbePairZWindow
        )
        self.events["nMuonVetoTagProbePairZWindowPass"] = ak.num(
            self.events.MuonVetoTagProbePairZWindowPass
        )
        self.events["nMuonVetoTagProbePairZWindowFail"] = ak.num(
            self.events.MuonVetoTagProbePairZWindowFail
        )
        self.events["nIsoTrackSearchPreMissingOuter"] = ak.num(
            self.events.IsoTrackSearchPreMissingOuter
        )
        self.events["nIsoTrackSearchPreLeptonVeto"] = ak.num(
            self.events.IsoTrackSearchPreLeptonVeto
        )
        self.events["nIsoTrackSearch"] = ak.num(self.events.IsoTrackSearch)

        track_diagnostics = {}
        diagnostics = {}
        for name, mask in search_track_cutflow_masks(self.events.IsoTrack).items():
            n_name = f"n{name[0].upper()}{name[1:]}"
            self.events[n_name] = ak.num(self.events.IsoTrack[mask])
            track_diagnostics[name] = self.events[n_name] >= 1

        diagnostics.update(track_diagnostics)

        event_diagnostics = search_event_cutflow_masks(self.events.AnalysisEvent)
        diagnostics.update(event_diagnostics)
        event_search_kinematics = event_diagnostics["event_dijetDphi2p5"]
        for name, mask in track_diagnostics.items():
            diagnostics[f"eventKinematics_{name}"] = event_search_kinematics & mask
        self.events["SearchDiag"] = ak.zip(diagnostics)

    def define_common_variables_before_presel(self, variation):
        pass
