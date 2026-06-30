"""PocketCoffea workflow for the custom disappearing-track NanoAOD."""

from __future__ import annotations

import awkward as ak

from pocket_coffea.workflows.base import BaseProcessorABC

from disapptrks.selections import (
    add_event_derived_fields,
    add_isotrack_derived_fields,
    base_probe_track_mask,
    search_event_cutflow_masks,
    search_track_cutflow_masks,
    search_track_mask,
)


class DisappTrksProcessor(BaseProcessorABC):
    def apply_object_preselection(self, variation):
        self.events["IsoTrack"] = add_isotrack_derived_fields(self.events)
        self.events["AnalysisEvent"] = add_event_derived_fields(self.events)
        self.events["IsoTrackProbe"] = self.events.IsoTrack[
            base_probe_track_mask(self.events.IsoTrack)
        ]
        search_diagnostic_masks = search_track_cutflow_masks(self.events.IsoTrack)
        self.events["IsoTrackSearchPreMissingOuter"] = self.events.IsoTrack[
            search_diagnostic_masks["track_calo10"]
        ]
        self.events["IsoTrackSearch"] = self.events.IsoTrack[
            search_track_mask(self.events.IsoTrack)
        ]

    def count_objects(self, variation):
        self.events["nIsoTrack"] = ak.num(self.events.IsoTrack)
        self.events["nIsoTrackProbe"] = ak.num(self.events.IsoTrackProbe)
        self.events["nIsoTrackSearchPreMissingOuter"] = ak.num(
            self.events.IsoTrackSearchPreMissingOuter
        )
        self.events["nIsoTrackSearch"] = ak.num(self.events.IsoTrackSearch)

        diagnostics = {}
        for name, mask in search_track_cutflow_masks(self.events.IsoTrack).items():
            n_name = f"n{name[0].upper()}{name[1:]}"
            self.events[n_name] = ak.num(self.events.IsoTrack[mask])
            diagnostics[name] = self.events[n_name] >= 1

        diagnostics.update(search_event_cutflow_masks(self.events.AnalysisEvent))
        self.events["SearchDiag"] = ak.zip(diagnostics)

    def define_common_variables_before_presel(self, variation):
        pass
