"""PocketCoffea workflow for the custom disappearing-track NanoAOD."""

from __future__ import annotations

import awkward as ak

from pocket_coffea.workflows.base import BaseProcessorABC

from disapptrks.selections import (
    add_event_derived_fields,
    add_isotrack_derived_fields,
    base_probe_track_mask,
    search_track_mask,
)


class DisappTrksProcessor(BaseProcessorABC):
    def apply_object_preselection(self, variation):
        self.events["IsoTrack"] = add_isotrack_derived_fields(self.events)
        self.events["AnalysisEvent"] = add_event_derived_fields(self.events)
        self.events["IsoTrackProbe"] = self.events.IsoTrack[
            base_probe_track_mask(self.events.IsoTrack)
        ]
        self.events["IsoTrackSearch"] = self.events.IsoTrack[
            search_track_mask(self.events.IsoTrack)
        ]

    def count_objects(self, variation):
        self.events["nIsoTrackProbe"] = ak.num(self.events.IsoTrackProbe)
        self.events["nIsoTrackSearch"] = ak.num(self.events.IsoTrackSearch)

    def define_common_variables_before_presel(self, variation):
        pass
