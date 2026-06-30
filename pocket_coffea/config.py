"""Minimal PocketCoffea configuration for local custom-NanoAOD validation.

Run from this directory after installing this package and the sibling
PocketCoffea checkout.
"""

from __future__ import annotations

import os

import cloudpickle

from pocket_coffea.parameters import defaults
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.lib.hist_manager import Axis, HistConf
from pocket_coffea.utils.configurator import Configurator

import cuts
import workflow
from cuts import has_disappearing_track, search_diagnostic_cuts, search_kinematics
from workflow import DisappTrksProcessor

cloudpickle.register_pickle_by_value(cuts)
cloudpickle.register_pickle_by_value(workflow)

localdir = os.path.dirname(os.path.abspath(__file__))
parameters = defaults.get_default_parameters()
diagnostic_categories = {
    f"diag_{name}": [cut] for name, cut in search_diagnostic_cuts.items()
}

cfg = Configurator(
    parameters=parameters,
    datasets={
        "jsons": [f"{localdir}/datasets/local_2024F_muon.json"],
        "filter": {"samples": ["DATA_Muon"], "year": ["2024"]},
    },
    workflow=DisappTrksProcessor,
    calibrators=[],
    skim=[],
    preselections=[],
    categories={
        "inclusive": [passthrough],
        "search": [search_kinematics, has_disappearing_track],
        **diagnostic_categories,
    },
    weights={"common": {"inclusive": []}, "bysample": {}},
    weights_classes=[],
    variations={"weights": {"common": {"inclusive": []}}},
    variables={
        "nIsoTrack": HistConf(
            [
                Axis(
                    coll="events",
                    field="nIsoTrack",
                    bins=20,
                    start=0,
                    stop=20,
                    label="N(isolated tracks)",
                )
            ]
        ),
        "nIsoTrackProbe": HistConf(
            [
                Axis(
                    coll="events",
                    field="nIsoTrackProbe",
                    bins=10,
                    start=0,
                    stop=10,
                    label="N(probe tracks)",
                )
            ]
        ),
        "nIsoTrackSearch": HistConf(
            [
                Axis(
                    coll="events",
                    field="nIsoTrackSearch",
                    bins=10,
                    start=0,
                    stop=10,
                    label="N(search tracks)",
                )
            ]
        ),
        "nIsoTrackSearchPreMissingOuter": HistConf(
            [
                Axis(
                    coll="events",
                    field="nIsoTrackSearchPreMissingOuter",
                    bins=10,
                    start=0,
                    stop=10,
                    label="N(search tracks before missing outer hit cut)",
                )
            ]
        ),
        "metNoMu_pt": HistConf(
            [
                Axis(
                    coll="AnalysisEvent",
                    field="METNoMu_pt",
                    bins=60,
                    start=0,
                    stop=600,
                    label=r"$p_T^{miss,no\ \mu}$ [GeV]",
                )
            ]
        ),
        "leadingJet_pt": HistConf(
            [
                Axis(
                    coll="AnalysisEvent",
                    field="leadingJet_pt",
                    bins=60,
                    start=0,
                    stop=600,
                    label="Leading selected jet $p_T$ [GeV]",
                )
            ]
        ),
        "leadingJet_metNoMu_deltaPhi": HistConf(
            [
                Axis(
                    coll="AnalysisEvent",
                    field="leadingJetMETNoMuDeltaPhi",
                    bins=32,
                    start=0,
                    stop=3.2,
                    label=r"$\Delta\phi$(leading jet, no-$\mu$ MET)",
                )
            ]
        ),
        "dijet_max_deltaPhi": HistConf(
            [
                Axis(
                    coll="AnalysisEvent",
                    field="dijetMaxDeltaPhi",
                    bins=34,
                    start=-0.2,
                    stop=3.2,
                    label=r"max $\Delta\phi$(selected dijets)",
                )
            ]
        ),
        "probeTrack_pt": HistConf(
            [
                Axis(
                    coll="IsoTrackProbe",
                    field="pt",
                    bins=60,
                    start=0,
                    stop=300,
                    label="Probe track $p_T$ [GeV]",
                )
            ]
        ),
        "probeTrack_eta": HistConf(
            [
                Axis(
                    coll="IsoTrackProbe",
                    field="eta",
                    bins=50,
                    start=-2.5,
                    stop=2.5,
                    label="Probe track eta",
                )
            ]
        ),
        "probeTrack_phi": HistConf(
            [
                Axis(
                    coll="IsoTrackProbe",
                    field="phi",
                    bins=64,
                    start=-3.2,
                    stop=3.2,
                    label="Probe track phi",
                )
            ]
        ),
        "probeTrack_caloEnergy": HistConf(
            [
                Axis(
                    coll="IsoTrackProbe",
                    field="caloEnergy",
                    bins=60,
                    start=0,
                    stop=60,
                    label="Probe track calorimeter energy [GeV]",
                )
            ]
        ),
        "probeTrack_dRMinJet": HistConf(
            [
                Axis(
                    coll="IsoTrackProbe",
                    field="dRMinJet",
                    bins=60,
                    start=-1,
                    stop=5,
                    label=r"Probe track min $\Delta R$(jet)",
                )
            ]
        ),
        "allTrack_missingOuterHits": HistConf(
            [
                Axis(
                    coll="IsoTrack",
                    field="missingOuterHits",
                    bins=21,
                    start=-0.5,
                    stop=20.5,
                    label="All isolated-track missing outer hits",
                )
            ]
        ),
        "probeTrack_missingOuterHits": HistConf(
            [
                Axis(
                    coll="IsoTrackProbe",
                    field="missingOuterHits",
                    bins=21,
                    start=-0.5,
                    stop=20.5,
                    label="Probe-track missing outer hits",
                )
            ]
        ),
        "preMissingOuterTrack_missingOuterHits": HistConf(
            [
                Axis(
                    coll="IsoTrackSearchPreMissingOuter",
                    field="missingOuterHits",
                    bins=21,
                    start=-0.5,
                    stop=20.5,
                    label="Search preselection track missing outer hits",
                )
            ]
        ),
        "searchTrack_pt": HistConf(
            [
                Axis(
                    coll="IsoTrackSearch",
                    field="pt",
                    bins=60,
                    start=0,
                    stop=300,
                    label="Search track $p_T$ [GeV]",
                )
            ]
        ),
        "searchTrack_eta_phi": HistConf(
            [
                Axis(
                    coll="IsoTrackSearch",
                    field="eta",
                    bins=50,
                    start=-2.5,
                    stop=2.5,
                    label="Search track eta",
                ),
                Axis(
                    coll="IsoTrackSearch",
                    field="phi",
                    bins=64,
                    start=-3.2,
                    stop=3.2,
                    label="Search track phi",
                ),
            ]
        ),
    },
    columns={},
)
