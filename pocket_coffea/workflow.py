"""PocketCoffea workflow for the custom disappearing-track NanoAOD."""

from __future__ import annotations

import awkward as ak

from pocket_coffea.workflows.base import BaseProcessorABC

from disapptrks.selections import (
    add_event_derived_fields,
    add_isotrack_derived_fields,
    add_muon_derived_fields,
    base_probe_track_mask,
    build_muon_veto_tag_probe_pairs,
    mass10_muon_probe_pair_mask,
    muon_tag_progression_masks,
    muon_tag_mask,
    muon_pveto_pair_pass_mask,
    muon_probe_pair_layer_mask,
    muon_veto_probe_track_mask,
    muon_veto_probe_track_cutflow_masks,
    muon_veto_pair_fail_mask,
    muon_veto_pair_pass_mask,
    os_mass10_muon_probe_pair_mask,
    os_muon_probe_pair_mask,
    os_z_window_muon_probe_pair_mask,
    search_event_cutflow_masks,
    search_track_cutflow_masks,
    search_track_mask,
    ss_mass10_muon_probe_pair_mask,
    ss_muon_probe_pair_mask,
    ss_z_window_muon_probe_pair_mask,
    z_window_muon_probe_pair_mask,
)

PVETO_LAYERS = ("NLayers4", "NLayers5", "NLayers6plus")


class DisappTrksProcessor(BaseProcessorABC):
    def apply_object_preselection(self, variation):
        self.events["Muon"] = add_muon_derived_fields(self.events)
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
        self.events["MuonVetoTagProbePairMass10"] = muon_veto_pairs[
            mass10_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairOSMass10"] = muon_veto_pairs[
            os_mass10_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairSS"] = muon_veto_pairs[
            ss_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairSSMass10"] = muon_veto_pairs[
            ss_mass10_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairZWindow"] = muon_veto_pairs[
            z_window_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairOSZWindow"] = muon_veto_pairs[
            os_z_window_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairZWindowPass"] = muon_veto_pairs[
            os_z_window_muon_probe_pair_mask(muon_veto_pairs)
            & muon_veto_pair_pass_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairZWindowFail"] = muon_veto_pairs[
            os_z_window_muon_probe_pair_mask(muon_veto_pairs)
            & muon_veto_pair_fail_mask(muon_veto_pairs)
        ]
        self.events["MuonPVetoTagProbePairZWindowPass"] = muon_veto_pairs[
            os_z_window_muon_probe_pair_mask(muon_veto_pairs)
            & muon_pveto_pair_pass_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairSSZWindow"] = muon_veto_pairs[
            ss_z_window_muon_probe_pair_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairSSZWindowPass"] = muon_veto_pairs[
            ss_z_window_muon_probe_pair_mask(muon_veto_pairs)
            & muon_veto_pair_pass_mask(muon_veto_pairs)
        ]
        self.events["MuonVetoTagProbePairSSZWindowFail"] = muon_veto_pairs[
            ss_z_window_muon_probe_pair_mask(muon_veto_pairs)
            & muon_veto_pair_fail_mask(muon_veto_pairs)
        ]
        self.events["MuonPVetoTagProbePairSSZWindowPass"] = muon_veto_pairs[
            ss_z_window_muon_probe_pair_mask(muon_veto_pairs)
            & muon_pveto_pair_pass_mask(muon_veto_pairs)
        ]
        for layer in PVETO_LAYERS:
            layer_mask = muon_probe_pair_layer_mask(muon_veto_pairs, layer)
            self.events[f"MuonVetoTagProbePairZWindow_{layer}"] = muon_veto_pairs[
                os_z_window_muon_probe_pair_mask(muon_veto_pairs) & layer_mask
            ]
            self.events[f"MuonPVetoTagProbePairZWindowPass_{layer}"] = muon_veto_pairs[
                os_z_window_muon_probe_pair_mask(muon_veto_pairs)
                & layer_mask
                & muon_pveto_pair_pass_mask(muon_veto_pairs)
            ]
            self.events[f"MuonVetoTagProbePairSSZWindow_{layer}"] = muon_veto_pairs[
                ss_z_window_muon_probe_pair_mask(muon_veto_pairs) & layer_mask
            ]
            self.events[f"MuonPVetoTagProbePairSSZWindowPass_{layer}"] = muon_veto_pairs[
                ss_z_window_muon_probe_pair_mask(muon_veto_pairs)
                & layer_mask
                & muon_pveto_pair_pass_mask(muon_veto_pairs)
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
        self.events["nMuonVetoTagProbePairMass10"] = ak.num(
            self.events.MuonVetoTagProbePairMass10
        )
        self.events["nMuonVetoTagProbePairOSMass10"] = ak.num(
            self.events.MuonVetoTagProbePairOSMass10
        )
        self.events["nMuonVetoTagProbePairSS"] = ak.num(
            self.events.MuonVetoTagProbePairSS
        )
        self.events["nMuonVetoTagProbePairSSMass10"] = ak.num(
            self.events.MuonVetoTagProbePairSSMass10
        )
        self.events["nMuonVetoTagProbePairZWindow"] = ak.num(
            self.events.MuonVetoTagProbePairZWindow
        )
        self.events["nMuonVetoTagProbePairOSZWindow"] = ak.num(
            self.events.MuonVetoTagProbePairOSZWindow
        )
        self.events["nMuonVetoTagProbePairZWindowPass"] = ak.num(
            self.events.MuonVetoTagProbePairZWindowPass
        )
        self.events["nMuonVetoTagProbePairZWindowFail"] = ak.num(
            self.events.MuonVetoTagProbePairZWindowFail
        )
        self.events["nMuonPVetoTagProbePairZWindowPass"] = ak.num(
            self.events.MuonPVetoTagProbePairZWindowPass
        )
        self.events["nMuonVetoTagProbePairSSZWindow"] = ak.num(
            self.events.MuonVetoTagProbePairSSZWindow
        )
        self.events["nMuonVetoTagProbePairSSZWindowPass"] = ak.num(
            self.events.MuonVetoTagProbePairSSZWindowPass
        )
        self.events["nMuonVetoTagProbePairSSZWindowFail"] = ak.num(
            self.events.MuonVetoTagProbePairSSZWindowFail
        )
        self.events["nMuonPVetoTagProbePairSSZWindowPass"] = ak.num(
            self.events.MuonPVetoTagProbePairSSZWindowPass
        )
        for layer in PVETO_LAYERS:
            self.events[f"nMuonVetoTagProbePairZWindow_{layer}"] = ak.num(
                self.events[f"MuonVetoTagProbePairZWindow_{layer}"]
            )
            self.events[f"nMuonPVetoTagProbePairZWindowPass_{layer}"] = ak.num(
                self.events[f"MuonPVetoTagProbePairZWindowPass_{layer}"]
            )
            self.events[f"nMuonVetoTagProbePairSSZWindow_{layer}"] = ak.num(
                self.events[f"MuonVetoTagProbePairSSZWindow_{layer}"]
            )
            self.events[f"nMuonPVetoTagProbePairSSZWindowPass_{layer}"] = ak.num(
                self.events[f"MuonPVetoTagProbePairSSZWindowPass_{layer}"]
            )
        self.events["nIsoTrackSearchPreMissingOuter"] = ak.num(
            self.events.IsoTrackSearchPreMissingOuter
        )
        self.events["nIsoTrackSearchPreLeptonVeto"] = ak.num(
            self.events.IsoTrackSearchPreLeptonVeto
        )
        self.events["nIsoTrackSearch"] = ak.num(self.events.IsoTrackSearch)

        muon_table16_diagnostics = {"event_singlemu_trigger": self.events.HLT.IsoMu24}
        muon_tag_masks = muon_tag_progression_masks(self.events.Muon)
        for name, mask in muon_tag_masks.items():
            self.events[f"n{name[0].upper()}{name[1:]}"] = ak.num(
                self.events.Muon[mask]
            )
            muon_table16_diagnostics[name] = self.events[
                f"n{name[0].upper()}{name[1:]}"
            ] >= 1

        has_selected_muon_tag = muon_table16_diagnostics["muon_selected_tag"]
        table16_track_masks = muon_veto_probe_track_cutflow_masks(self.events.IsoTrack)
        for name, mask in table16_track_masks.items():
            self.events[f"n{name[0].upper()}{name[1:]}Table16"] = ak.num(
                self.events.IsoTrack[mask]
            )
            muon_table16_diagnostics[name] = (
                has_selected_muon_tag
                & (self.events[f"n{name[0].upper()}{name[1:]}Table16"] >= 1)
            )

        table16_mass_probe_tracks = self.events.IsoTrack[
            table16_track_masks["track_dRJet0p5"]
        ]
        table16_mass_pairs = build_muon_veto_tag_probe_pairs(
            self.events.MuonTag, table16_mass_probe_tracks
        )
        table16_probe_tracks = self.events.IsoTrack[table16_track_masks["track_calo10"]]
        table16_pairs = build_muon_veto_tag_probe_pairs(
            self.events.MuonTag, table16_probe_tracks
        )
        muon_table16_diagnostics.update(
            {
                "pair_mass10": ak.num(
                    table16_mass_pairs[
                        mass10_muon_probe_pair_mask(table16_mass_pairs)
                    ]
                )
                >= 1,
                "pair_zwindow": ak.num(
                    table16_pairs[z_window_muon_probe_pair_mask(table16_pairs)]
                )
                >= 1,
                "pair_os": ak.num(
                    table16_pairs[os_z_window_muon_probe_pair_mask(table16_pairs)]
                )
                >= 1,
            }
        )
        for layer in PVETO_LAYERS:
            layer_mask = muon_probe_pair_layer_mask(table16_pairs, layer)
            muon_table16_diagnostics[f"layer_{layer}"] = ak.num(
                table16_pairs[os_z_window_muon_probe_pair_mask(table16_pairs) & layer_mask]
            ) >= 1
        self.events["MuonTable16Diag"] = ak.zip(muon_table16_diagnostics)

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
