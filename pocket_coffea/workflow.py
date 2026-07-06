"""PocketCoffea workflow for the custom disappearing-track NanoAOD."""

from __future__ import annotations

import os
from pathlib import Path

import awkward as ak
import numpy as np

from pocket_coffea.workflows.base import BaseProcessorABC

from disapptrks.selections import (
    add_event_derived_fields,
    add_isotrack_derived_fields,
    add_muon_derived_fields,
    base_probe_track_mask,
    build_lepton_veto_tag_probe_pairs,
    build_muon_veto_tag_probe_pairs,
    electron_pveto_pair_pass_mask,
    electron_tag_progression_masks,
    electron_tag_mask,
    fake_track_no_d0_mask,
    generic_probe_pair_layer_mask,
    lepton_veto_probe_track_mask,
    low_mt_mask,
    mass10_muon_probe_pair_mask,
    mass_window_pair_mask,
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
    os_mass_window_pair_mask,
    os_z_window_muon_probe_pair_mask,
    search_event_cutflow_masks,
    search_track_cutflow_masks,
    search_track_mask,
    single_electron_trigger_mask,
    ss_mass10_muon_probe_pair_mask,
    ss_muon_probe_pair_mask,
    ss_mass_window_pair_mask,
    ss_z_window_muon_probe_pair_mask,
    tau_pveto_pair_pass_mask,
    z_window_muon_probe_pair_mask,
)

PVETO_LAYERS = ("NLayers4", "NLayers5", "NLayers6plus")
ELECTRON_MASS = 0.000511
MUON_MASS = 0.105658

JET_VETO_MAP_FILES = {
    "2022_preEE": "Run3-22CDSep23-Summer22-NanoAODv12_jetvetomaps.json.gz",
    "2022_postEE": "Run3-22EFGSep23-Summer22EE-NanoAODv12_jetvetomaps.json.gz",
    "2023_preBPix": "Run3-23CSep23-Summer23-NanoAODv12_jetvetomaps.json.gz",
    "2023_postBPix": "Run3-23DSep23-Summer23BPix-NanoAODv12_jetvetomaps.json.gz",
    "2024": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15_jetvetomaps.json.gz",
    "2025": "Run3-25Prompt-Winter25-NanoAODv15_jetvetomaps.json.gz",
}


def _all_true_like(events):
    return ak.ones_like(events.HLT.IsoMu24, dtype=bool)


def _event_flag(events, name: str, *, default: bool = True):
    if "Flag" in events.fields and name in events.Flag.fields:
        return events.Flag[name]
    if name in events.fields:
        return events[name]
    return _all_true_like(events) if default else ~_all_true_like(events)


def _met_filters_mask(events):
    if "Flag" in events.fields and "METFilters" in events.Flag.fields:
        return events.Flag.METFilters
    if "METFilters" in events.fields:
        return events.METFilters

    filters = [
        "goodVertices",
        "globalSuperTightHalo2016Filter",
        "EcalDeadCellTriggerPrimitiveFilter",
        "BadPFMuonFilter",
        "BadPFMuonDzFilter",
        "hfNoisyHitsFilter",
        "eeBadScFilter",
        "ecalBadCalibFilter",
    ]
    mask = _all_true_like(events)
    for name in filters:
        mask = mask & _event_flag(events, name, default=True)
    return mask


def _met_for_transverse_mass(events):
    """Return the event-level MET collection used for lepton MT cuts.

    Some of the newer/custom NanoAOD files used in the 2024/2025 productions do
    not store the classic ``MET`` collection, but do store ``PuppiMET``.  The
    low-MT tag selection only needs an event-level ``pt`` and ``phi``, so use
    ``MET`` when available and fall back to ``PuppiMET`` otherwise.
    """
    for name in ("MET", "PuppiMET"):
        if name in events.fields:
            return events[name]
    raise AttributeError(
        "No MET-like collection found for lepton transverse-mass cuts; "
        "expected either MET or PuppiMET."
    )


def _jet_veto_map_parameter_year(year, era, processor_params):
    year = str(year)
    for container in (
        processor_params.jet_scale_factors.vetomaps,
        processor_params.lumi.goldenJSON,
        processor_params.event_flags,
    ):
        if year in container:
            return year

    if year == "2022":
        return "2022_preEE" if era in ("C", "D") else "2022_postEE"
    if year == "2023":
        return "2023_preBPix" if era == "C" else "2023_postBPix"
    return year


def _local_jet_veto_map_path(mapped_year):
    filename = JET_VETO_MAP_FILES.get(str(mapped_year))
    if filename is None:
        return None

    search_dirs = []
    env_dir = os.environ.get("DISAPPTRKS_JET_VETO_MAP_DIR")
    if env_dir:
        search_dirs.append(Path(env_dir))
    search_dirs.append(Path(__file__).resolve().parent / "data" / "jet_veto_maps")
    search_dirs.append(Path.cwd() / "data" / "jet_veto_maps")
    search_dirs.append(Path.cwd() / "jet_veto_maps")

    for directory in search_dirs:
        candidates = [
            directory / filename,
            directory / str(mapped_year) / "jetvetomaps.json.gz",
            directory / filename.removesuffix("_jetvetomaps.json.gz") / "jetvetomaps.json.gz",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
    return None


def _evaluate_jet_veto_map(events, processor_params, mapped_year, payload_file):
    import correctionlib
    from pocket_coffea.lib.jets import compute_jetId

    jets = ak.with_field(
        events["Jet"],
        compute_jetId(events, "Jet", processor_params, mapped_year),
        "jetId_corrected",
    )
    mask_for_veto_map = (
        (jets["jetId_corrected"] >= 6)
        & (abs(jets.eta) < 5.19)
        & (jets.pt > 15.0)
        & ((jets["neEmEF"] + jets["chEmEF"]) < 0.9)
    )
    jets = jets[mask_for_veto_map]

    cset = correctionlib.CorrectionSet.from_file(str(payload_file))
    corr = cset[processor_params.jet_scale_factors.vetomaps[mapped_year]["name"]]
    eta_flat = ak.to_numpy(ak.flatten(jets.eta))
    phi_flat = np.clip(ak.to_numpy(ak.flatten(jets.phi)), -3.14159, 3.14159)
    eta_counts = ak.num(jets.eta)
    weight = ak.unflatten(
        corr.evaluate("jetvetomap", eta_flat, phi_flat),
        counts=eta_counts,
    )
    event_mask = ak.sum(weight, axis=-1) == 0
    return ak.where(ak.is_none(event_mask), False, event_mask)


def _golden_json_mask(
    events,
    *,
    processor_params=None,
    year=None,
    era=None,
    sample=None,
    is_mc=False,
):
    if is_mc or processor_params is None or year is None:
        return _all_true_like(events)

    if "run" not in events.fields or "luminosityBlock" not in events.fields:
        return _all_true_like(events)

    from pocket_coffea.lib.cut_functions import apply_golden_json

    golden_json_year = _jet_veto_map_parameter_year(year, era, processor_params)
    return apply_golden_json(
        events,
        params={},
        year=golden_json_year,
        processor_params=processor_params,
        sample=sample,
        isMC=is_mc,
    )


def _jet_veto_map_mask(
    events,
    *,
    processor_params=None,
    year=None,
    era=None,
    sample=None,
    is_mc=False,
):
    # The JME jet-veto-map decision is not guaranteed to exist in central
    # NanoAOD.  Legacy DisappTrks computed and saved it as ``jetVeto2022``; some
    # custom NanoAOD productions may store the same decision under Flag.  When
    # it is absent, use PocketCoffea's correctionlib-based Run-3 jet-veto-map
    # evaluator, preferring local payloads in environments without CMS CVMFS.
    if "Flag" in events.fields and "jetVeto2022" in events.Flag.fields:
        return events.Flag.jetVeto2022
    if "jetVeto2022" in events.fields:
        return events.jetVeto2022

    if processor_params is None or year is None:
        if os.environ.get("DISAPPTRKS_ALLOW_MISSING_JET_VETO_MAP", "").lower() in (
            "1",
            "true",
            "yes",
            "on",
        ):
            return _all_true_like(events)
        raise RuntimeError("Cannot apply the JME jet-veto map without processor parameters and year.")

    veto_map_year = _jet_veto_map_parameter_year(year, era, processor_params)
    local_payload = _local_jet_veto_map_path(veto_map_year)
    if local_payload is not None:
        return _evaluate_jet_veto_map(events, processor_params, veto_map_year, local_payload)

    try:
        from pocket_coffea.lib.cut_functions import get_JetVetoMap_Mask

        return get_JetVetoMap_Mask(
            events,
            params={},
            year=veto_map_year,
            processor_params=processor_params,
            sample=sample,
            isMC=is_mc,
        )
    except Exception as exc:
        if os.environ.get("DISAPPTRKS_ALLOW_MISSING_JET_VETO_MAP", "").lower() in (
            "1",
            "true",
            "yes",
            "on",
        ):
            return _all_true_like(events)
        raise RuntimeError(
            "Could not apply the JME jet-veto map. Provide local payloads in "
            "pocket_coffea/data/jet_veto_maps, set DISAPPTRKS_JET_VETO_MAP_DIR, "
            "or explicitly set DISAPPTRKS_ALLOW_MISSING_JET_VETO_MAP=1 for a "
            "non-production diagnostic run."
        ) from exc


class DisappTrksProcessor(BaseProcessorABC):
    def _store_lepton_pveto_pairs(
        self,
        *,
        prefix,
        tags,
        probes,
        tag_mass,
        probe_mass,
        window_low,
        window_high,
        pass_mask_function,
    ):
        pairs = build_lepton_veto_tag_probe_pairs(
            tags,
            probes,
            tag_mass=tag_mass,
            probe_mass=probe_mass,
        )
        self.events[f"{prefix}TagProbePair"] = pairs
        self.events[f"{prefix}TagProbePairMassWindow"] = pairs[
            mass_window_pair_mask(pairs, window_low, window_high)
        ]
        self.events[f"{prefix}TagProbePairOSMassWindow"] = pairs[
            os_mass_window_pair_mask(pairs, window_low, window_high)
        ]
        self.events[f"{prefix}TagProbePairSSMassWindow"] = pairs[
            ss_mass_window_pair_mask(pairs, window_low, window_high)
        ]
        self.events[f"{prefix}PVetoTagProbePairMassWindowPass"] = pairs[
            os_mass_window_pair_mask(pairs, window_low, window_high)
            & pass_mask_function(pairs)
        ]
        self.events[f"{prefix}PVetoTagProbePairSSMassWindowPass"] = pairs[
            ss_mass_window_pair_mask(pairs, window_low, window_high)
            & pass_mask_function(pairs)
        ]
        for layer in PVETO_LAYERS:
            layer_mask = generic_probe_pair_layer_mask(pairs, layer)
            self.events[f"{prefix}TagProbePairMassWindow_{layer}"] = pairs[
                os_mass_window_pair_mask(pairs, window_low, window_high) & layer_mask
            ]
            self.events[f"{prefix}PVetoTagProbePairMassWindowPass_{layer}"] = pairs[
                os_mass_window_pair_mask(pairs, window_low, window_high)
                & layer_mask
                & pass_mask_function(pairs)
            ]
            self.events[f"{prefix}TagProbePairSSMassWindow_{layer}"] = pairs[
                ss_mass_window_pair_mask(pairs, window_low, window_high) & layer_mask
            ]
            self.events[f"{prefix}PVetoTagProbePairSSMassWindowPass_{layer}"] = pairs[
                ss_mass_window_pair_mask(pairs, window_low, window_high)
                & layer_mask
                & pass_mask_function(pairs)
            ]

    def apply_object_preselection(self, variation):
        self.events["Muon"] = add_muon_derived_fields(self.events)
        self.events["IsoTrack"] = add_isotrack_derived_fields(self.events)
        self.events["AnalysisEvent"] = add_event_derived_fields(self.events)
        tag_met = _met_for_transverse_mass(self.events)
        self.events["MuonTag"] = self.events.Muon[muon_tag_mask(self.events.Muon)]
        self.events["MuonLowMTTag"] = self.events.MuonTag[
            low_mt_mask(self.events.MuonTag, tag_met)
        ]
        self.events["ElectronTag"] = self.events.Electron[
            electron_tag_mask(self.events.Electron, self.events)
        ]
        self.events["ElectronLowMTTag"] = self.events.ElectronTag[
            low_mt_mask(self.events.ElectronTag, tag_met)
        ]
        self.events["IsoTrackProbe"] = self.events.IsoTrack[
            base_probe_track_mask(self.events.IsoTrack)
        ]
        self.events["MuonVetoProbeTrack"] = self.events.IsoTrack[
            muon_veto_probe_track_mask(self.events.IsoTrack)
        ]
        self.events["ElectronVetoProbeTrack"] = self.events.IsoTrack[
            lepton_veto_probe_track_mask(self.events.IsoTrack, measured_veto="electron")
        ]
        self.events["TauVetoProbeTrack"] = self.events.IsoTrack[
            lepton_veto_probe_track_mask(self.events.IsoTrack, measured_veto="tau")
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
        self._store_lepton_pveto_pairs(
            prefix="Electron",
            tags=self.events.ElectronTag,
            probes=self.events.ElectronVetoProbeTrack,
            tag_mass=ELECTRON_MASS,
            probe_mass=ELECTRON_MASS,
            window_low=91.1876 - 10.0,
            window_high=91.1876 + 10.0,
            pass_mask_function=electron_pveto_pair_pass_mask,
        )
        self._store_lepton_pveto_pairs(
            prefix="TauMu",
            tags=self.events.MuonLowMTTag,
            probes=self.events.TauVetoProbeTrack,
            tag_mass=MUON_MASS,
            probe_mass=MUON_MASS,
            window_low=91.1876 - 50.0,
            window_high=91.1876 - 15.0,
            pass_mask_function=tau_pveto_pair_pass_mask,
        )
        self._store_lepton_pveto_pairs(
            prefix="TauEle",
            tags=self.events.ElectronLowMTTag,
            probes=self.events.TauVetoProbeTrack,
            tag_mass=ELECTRON_MASS,
            probe_mass=ELECTRON_MASS,
            window_low=91.1876 - 50.0,
            window_high=91.1876 - 15.0,
            pass_mask_function=tau_pveto_pair_pass_mask,
        )
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
        self.events["nElectronTag"] = ak.num(self.events.ElectronTag)
        self.events["nMuonLowMTTag"] = ak.num(self.events.MuonLowMTTag)
        self.events["nElectronLowMTTag"] = ak.num(self.events.ElectronLowMTTag)
        self.events["nMuonVetoProbeTrack"] = ak.num(self.events.MuonVetoProbeTrack)
        self.events["nElectronVetoProbeTrack"] = ak.num(self.events.ElectronVetoProbeTrack)
        self.events["nTauVetoProbeTrack"] = ak.num(self.events.TauVetoProbeTrack)
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
        for prefix in ("Electron", "TauMu", "TauEle"):
            self.events[f"n{prefix}TagProbePair"] = ak.num(
                self.events[f"{prefix}TagProbePair"]
            )
            self.events[f"n{prefix}TagProbePairMassWindow"] = ak.num(
                self.events[f"{prefix}TagProbePairMassWindow"]
            )
            self.events[f"n{prefix}TagProbePairOSMassWindow"] = ak.num(
                self.events[f"{prefix}TagProbePairOSMassWindow"]
            )
            self.events[f"n{prefix}TagProbePairSSMassWindow"] = ak.num(
                self.events[f"{prefix}TagProbePairSSMassWindow"]
            )
            self.events[f"n{prefix}PVetoTagProbePairMassWindowPass"] = ak.num(
                self.events[f"{prefix}PVetoTagProbePairMassWindowPass"]
            )
            self.events[f"n{prefix}PVetoTagProbePairSSMassWindowPass"] = ak.num(
                self.events[f"{prefix}PVetoTagProbePairSSMassWindowPass"]
            )
            for layer in PVETO_LAYERS:
                self.events[f"n{prefix}TagProbePairMassWindow_{layer}"] = ak.num(
                    self.events[f"{prefix}TagProbePairMassWindow_{layer}"]
                )
                self.events[f"n{prefix}PVetoTagProbePairMassWindowPass_{layer}"] = ak.num(
                    self.events[f"{prefix}PVetoTagProbePairMassWindowPass_{layer}"]
                )
                self.events[f"n{prefix}TagProbePairSSMassWindow_{layer}"] = ak.num(
                    self.events[f"{prefix}TagProbePairSSMassWindow_{layer}"]
                )
                self.events[f"n{prefix}PVetoTagProbePairSSMassWindowPass_{layer}"] = ak.num(
                    self.events[f"{prefix}PVetoTagProbePairSSMassWindowPass_{layer}"]
                )
        self.events["nIsoTrackSearchPreMissingOuter"] = ak.num(
            self.events.IsoTrackSearchPreMissingOuter
        )
        self.events["nIsoTrackSearchPreLeptonVeto"] = ak.num(
            self.events.IsoTrackSearchPreLeptonVeto
        )
        self.events["nIsoTrackSearch"] = ak.num(self.events.IsoTrackSearch)

        fake_basic3hits_d0_signal = self.events.IsoTrack[
            fake_track_no_d0_mask(
                self.events.IsoTrack,
                layer="NLayers4",
                d0_region="signal",
            )
        ]
        fake_basic3hits_d0_sideband = self.events.IsoTrack[
            fake_track_no_d0_mask(
                self.events.IsoTrack,
                layer="NLayers4",
                d0_region="sideband",
            )
        ]
        self.events["nFakeBasic3HitsD0Signal"] = ak.num(fake_basic3hits_d0_signal)
        self.events["nFakeBasic3HitsD0Sideband"] = ak.num(fake_basic3hits_d0_sideband)
        for layer in (*PVETO_LAYERS, "combinedBins"):
            self.events[f"nFakeControl_{layer}"] = ak.num(
                self.events.IsoTrack[
                    fake_track_no_d0_mask(
                        self.events.IsoTrack,
                        layer=layer,
                        d0_region="sideband",
                    )
                ]
            )

        event_golden_json = _golden_json_mask(
            self.events,
            processor_params=self.params,
            year=self._year,
            era=self._era,
            sample=self._sample,
            is_mc=self._isMC,
        )
        event_singlemu_trigger = event_golden_json & self.events.HLT.IsoMu24
        event_met_filters = event_singlemu_trigger & _met_filters_mask(self.events)
        event_jet_veto_map = event_met_filters & _jet_veto_map_mask(
            self.events,
            processor_params=self.params,
            year=self._year,
            era=self._era,
            sample=self._sample,
            is_mc=self._isMC,
        )
        muon_table16_diagnostics = {
            "event_singlemu_trigger": event_singlemu_trigger,
            "event_met_filters": event_met_filters,
            "event_jet_veto_map": event_jet_veto_map,
        }
        muon_tag_masks = muon_tag_progression_masks(self.events.Muon)
        for name, mask in muon_tag_masks.items():
            self.events[f"n{name[0].upper()}{name[1:]}"] = ak.num(
                self.events.Muon[mask]
            )
            muon_table16_diagnostics[name] = (
                event_jet_veto_map
                & (self.events[f"n{name[0].upper()}{name[1:]}"] >= 1)
            )

        has_selected_muon_tag = muon_table16_diagnostics["muon_selected_tag"]
        table16_track_masks = muon_veto_probe_track_cutflow_masks(self.events.IsoTrack)
        pre_pair_track_fields = {
            "track_pt30",
            "track_eta2p1",
            "track_noDTWheelGap",
            "track_noECALCrack",
            "track_noCSCTransition",
            "track_fiducialECAL",
            "track_dzOrLambda",
            "track_pixelHits4",
            "track_noMissingInner",
            "track_noMissingMiddle",
            "track_chargedIso0p05",
            "track_dxy0p02",
            "track_dz0p5",
            "track_dRJet0p5",
        }
        for name, mask in table16_track_masks.items():
            self.events[f"n{name[0].upper()}{name[1:]}Table16"] = ak.num(
                self.events.IsoTrack[mask]
            )
            if name in pre_pair_track_fields:
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
        table16_electron_probe_tracks = self.events.IsoTrack[
            table16_track_masks["track_electronVeto"]
        ]
        table16_electron_pairs = build_muon_veto_tag_probe_pairs(
            self.events.MuonTag, table16_electron_probe_tracks
        )
        table16_tau_probe_tracks = self.events.IsoTrack[
            table16_track_masks["track_tauVeto"]
        ]
        table16_tau_pairs = build_muon_veto_tag_probe_pairs(
            self.events.MuonTag, table16_tau_probe_tracks
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
                "track_electronVeto": ak.num(
                    table16_electron_pairs[
                        mass10_muon_probe_pair_mask(table16_electron_pairs)
                    ]
                )
                >= 1,
                "track_tauVeto": ak.num(
                    table16_tau_pairs[
                        mass10_muon_probe_pair_mask(table16_tau_pairs)
                    ]
                )
                >= 1,
                "track_calo10": ak.num(
                    table16_pairs[
                        mass10_muon_probe_pair_mask(table16_pairs)
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
        muon_table16_diagnostics["track_probe_before_layer"] = (
            muon_table16_diagnostics["track_calo10"]
        )
        muon_table16_diagnostics["layer_combinedBins"] = ak.num(
            table16_pairs[
                os_z_window_muon_probe_pair_mask(table16_pairs)
                & muon_probe_pair_layer_mask(table16_pairs, "combinedBins")
            ]
        ) >= 1
        self.events["MuonTable16Diag"] = ak.zip(muon_table16_diagnostics)

        event_singleele_trigger = event_golden_json & single_electron_trigger_mask(
            self.events
        )
        event_ele_met_filters = event_singleele_trigger & _met_filters_mask(
            self.events
        )
        event_ele_jet_veto_map = event_ele_met_filters & _jet_veto_map_mask(
            self.events,
            processor_params=self.params,
            year=self._year,
            era=self._era,
            sample=self._sample,
            is_mc=self._isMC,
        )
        electron_pveto_diagnostics = {
            "event_singleele_trigger": event_singleele_trigger,
            "event_met_filters": event_ele_met_filters,
            "event_jet_veto_map": event_ele_jet_veto_map,
        }
        electron_tag_masks = electron_tag_progression_masks(
            self.events.Electron, self.events
        )
        for name, mask in electron_tag_masks.items():
            self.events[f"n{name[0].upper()}{name[1:]}"] = ak.num(
                self.events.Electron[mask]
            )
            electron_pveto_diagnostics[name] = (
                event_ele_jet_veto_map
                & (self.events[f"n{name[0].upper()}{name[1:]}"] >= 1)
            )

        has_selected_electron_tag = electron_pveto_diagnostics[
            "electron_selected_tag"
        ]
        electron_track_masks = muon_veto_probe_track_cutflow_masks(
            self.events.IsoTrack
        )
        for name in (
            "track_pt30",
            "track_eta2p1",
            "track_noDTWheelGap",
            "track_noECALCrack",
            "track_noCSCTransition",
            "track_fiducialECAL",
            "track_dzOrLambda",
            "track_pixelHits4",
            "track_noMissingInner",
            "track_noMissingMiddle",
            "track_chargedIso0p05",
            "track_dxy0p02",
            "track_dz0p5",
            "track_dRJet0p5",
        ):
            n_name = f"n{name[0].upper()}{name[1:]}ElectronPVeto"
            self.events[n_name] = ak.num(self.events.IsoTrack[electron_track_masks[name]])
            electron_pveto_diagnostics[name] = (
                has_selected_electron_tag & (self.events[n_name] >= 1)
            )

        electron_mass_probe_tracks = self.events.IsoTrack[
            electron_track_masks["track_dRJet0p5"]
        ]
        electron_mass_pairs = build_lepton_veto_tag_probe_pairs(
            self.events.ElectronTag,
            electron_mass_probe_tracks,
            tag_mass=ELECTRON_MASS,
            probe_mass=ELECTRON_MASS,
        )
        electron_muon_veto_mask = electron_track_masks["track_dRJet0p5"] & (
            (self.events.IsoTrack.dRMinMuon < 0.0)
            | (self.events.IsoTrack.dRMinMuon > 0.15)
        )
        electron_tau_veto_mask = electron_muon_veto_mask & (
            (self.events.IsoTrack.dRMinTauHad < 0.0)
            | (self.events.IsoTrack.dRMinTauHad > 0.15)
        )
        electron_calo_mask = electron_tau_veto_mask & (
            self.events.IsoTrack.caloEnergy < 10.0
        )
        electron_muon_veto_pairs = build_lepton_veto_tag_probe_pairs(
            self.events.ElectronTag,
            self.events.IsoTrack[electron_muon_veto_mask],
            tag_mass=ELECTRON_MASS,
            probe_mass=ELECTRON_MASS,
        )
        electron_tau_veto_pairs = build_lepton_veto_tag_probe_pairs(
            self.events.ElectronTag,
            self.events.IsoTrack[electron_tau_veto_mask],
            tag_mass=ELECTRON_MASS,
            probe_mass=ELECTRON_MASS,
        )
        electron_pairs = build_lepton_veto_tag_probe_pairs(
            self.events.ElectronTag,
            self.events.IsoTrack[electron_tau_veto_mask],
            tag_mass=ELECTRON_MASS,
            probe_mass=ELECTRON_MASS,
        )
        electron_calo_pairs = build_lepton_veto_tag_probe_pairs(
            self.events.ElectronTag,
            self.events.IsoTrack[electron_calo_mask],
            tag_mass=ELECTRON_MASS,
            probe_mass=ELECTRON_MASS,
        )
        electron_z_window = mass_window_pair_mask(
            electron_pairs, 91.1876 - 10.0, 91.1876 + 10.0
        )
        electron_os_z_window = os_mass_window_pair_mask(
            electron_pairs, 91.1876 - 10.0, 91.1876 + 10.0
        )
        electron_pveto_diagnostics.update(
            {
                "pair_mass10": ak.num(electron_mass_pairs[electron_mass_pairs.mass > 10.0])
                >= 1,
                "track_muonVeto": ak.num(
                    electron_muon_veto_pairs[electron_muon_veto_pairs.mass > 10.0]
                )
                >= 1,
                "track_tauVeto": ak.num(
                    electron_tau_veto_pairs[electron_tau_veto_pairs.mass > 10.0]
                )
                >= 1,
                "track_calo10": ak.num(
                    electron_calo_pairs[electron_calo_pairs.mass > 10.0]
                )
                >= 1,
                "track_probe_before_layer": ak.num(
                    electron_pairs[electron_pairs.mass > 10.0]
                )
                >= 1,
                "pair_zwindow": ak.num(electron_pairs[electron_z_window]) >= 1,
                "pair_os": ak.num(electron_pairs[electron_os_z_window]) >= 1,
                "layer_combinedBins": ak.num(
                    electron_pairs[
                        electron_os_z_window
                        & generic_probe_pair_layer_mask(electron_pairs, "combinedBins")
                    ]
                )
                >= 1,
                "pair_pass_electron_pveto": ak.num(
                    electron_pairs[
                        electron_os_z_window
                        & generic_probe_pair_layer_mask(electron_pairs, "combinedBins")
                        & electron_pveto_pair_pass_mask(electron_pairs)
                    ]
                )
                >= 1,
            }
        )
        self.events["ElectronPVetoDiag"] = ak.zip(electron_pveto_diagnostics)

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
