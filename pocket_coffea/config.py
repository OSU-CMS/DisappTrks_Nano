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
from cuts import (
    has_muon_tag,
    has_muon_veto_probe_track,
    has_muon_veto_tag_probe_pair,
    has_muon_veto_zwindow_pair,
    has_muon_veto_zwindow_pass_pair,
)
from workflow import DisappTrksProcessor

cloudpickle.register_pickle_by_value(cuts)
cloudpickle.register_pickle_by_value(workflow)

localdir = os.path.dirname(os.path.abspath(__file__))
parameters = defaults.get_default_parameters()
diagnostic_categories = {
    f"diag_{name}": [cut] for name, cut in search_diagnostic_cuts.items()
}
outer_hit_variants = {
    "missingOuterHits": "tracker layers without measurement",
    "hp_nLostHitsOuter": "lost hits",
    "hp_nLostTrackerHitsOuter": "lost tracker hits",
    "hp_nLostPixelHitsOuter": "lost pixel hits",
    "hp_nLostStripHitsOuter": "lost strip hits",
    "hp_pixelLayersWithoutMeasurementOuter": "pixel layers without measurement",
    "hp_stripLayersWithoutMeasurementOuter": "strip layers without measurement",
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
        "muon_veto_tag": [has_muon_tag],
        "muon_veto_probe": [has_muon_tag, has_muon_veto_probe_track],
        "muon_veto_pair": [has_muon_veto_tag_probe_pair],
        "muon_veto_zwindow": [has_muon_veto_zwindow_pair],
        "muon_veto_zwindow_pass": [has_muon_veto_zwindow_pass_pair],
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
        "nMuonTag": HistConf(
            [
                Axis(
                    coll="events",
                    field="nMuonTag",
                    bins=10,
                    start=0,
                    stop=10,
                    label="N(tight isolated muon tags)",
                )
            ]
        ),
        "nMuonVetoProbeTrack": HistConf(
            [
                Axis(
                    coll="events",
                    field="nMuonVetoProbeTrack",
                    bins=10,
                    start=0,
                    stop=10,
                    label="N(muon-veto probe tracks)",
                )
            ]
        ),
        "nMuonVetoTagProbePair": HistConf(
            [
                Axis(
                    coll="events",
                    field="nMuonVetoTagProbePair",
                    bins=20,
                    start=0,
                    stop=20,
                    label="N(muon tag-probe pairs)",
                )
            ]
        ),
        "nMuonVetoTagProbePairZWindow": HistConf(
            [
                Axis(
                    coll="events",
                    field="nMuonVetoTagProbePairZWindow",
                    bins=10,
                    start=0,
                    stop=10,
                    label="N(OS muon tag-probe pairs in Z window)",
                )
            ]
        ),
        "nMuonVetoTagProbePairZWindowPass": HistConf(
            [
                Axis(
                    coll="events",
                    field="nMuonVetoTagProbePairZWindowPass",
                    bins=10,
                    start=0,
                    stop=10,
                    label="N(OS Z-window pairs passing muon veto)",
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
        "nIsoTrackSearchPreLeptonVeto": HistConf(
            [
                Axis(
                    coll="events",
                    field="nIsoTrackSearchPreLeptonVeto",
                    bins=10,
                    start=0,
                    stop=10,
                    label="N(search tracks before lepton vetoes)",
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
        "muonTag_pt": HistConf(
            [
                Axis(
                    coll="MuonTag",
                    field="pt",
                    bins=60,
                    start=0,
                    stop=300,
                    label="Muon tag $p_T$ [GeV]",
                )
            ]
        ),
        "muonTag_eta": HistConf(
            [
                Axis(
                    coll="MuonTag",
                    field="eta",
                    bins=50,
                    start=-2.5,
                    stop=2.5,
                    label="Muon tag eta",
                )
            ]
        ),
        "muonVetoProbeTrack_pt": HistConf(
            [
                Axis(
                    coll="MuonVetoProbeTrack",
                    field="pt",
                    bins=60,
                    start=0,
                    stop=300,
                    label="Muon-veto probe track $p_T$ [GeV]",
                )
            ]
        ),
        "muonVetoProbeTrack_eta": HistConf(
            [
                Axis(
                    coll="MuonVetoProbeTrack",
                    field="eta",
                    bins=50,
                    start=-2.5,
                    stop=2.5,
                    label="Muon-veto probe track eta",
                )
            ]
        ),
        "muonVetoProbeTrack_dRMinMuon": HistConf(
            [
                Axis(
                    coll="MuonVetoProbeTrack",
                    field="dRMinMuon",
                    bins=60,
                    start=-1,
                    stop=5,
                    label=r"Probe track min $\Delta R$(muon)",
                )
            ]
        ),
        "muonVetoTagProbePair_mass": HistConf(
            [
                Axis(
                    coll="MuonVetoTagProbePair",
                    field="mass",
                    bins=80,
                    start=0,
                    stop=200,
                    label="Muon tag-probe mass [GeV]",
                )
            ]
        ),
        "muonVetoTagProbePairZWindow_mass": HistConf(
            [
                Axis(
                    coll="MuonVetoTagProbePairZWindow",
                    field="mass",
                    bins=40,
                    start=70,
                    stop=110,
                    label="OS Z-window muon tag-probe mass [GeV]",
                )
            ]
        ),
        "muonVetoTagProbePairZWindow_probe_dRMinMuon": HistConf(
            [
                Axis(
                    coll="MuonVetoTagProbePairZWindow",
                    field="probe_dRMinMuon",
                    bins=60,
                    start=-1,
                    stop=5,
                    label=r"OS Z-window probe min $\Delta R$(muon)",
                )
            ]
        ),
        **{
            f"allTrack_{field}": HistConf(
                [
                    Axis(
                        coll="IsoTrack",
                        field=field,
                        bins=21,
                        start=-0.5,
                        stop=20.5,
                        label=f"All isolated-track outer {label}",
                    )
                ]
            )
            for field, label in outer_hit_variants.items()
        },
        **{
            f"probeTrack_{field}": HistConf(
                [
                    Axis(
                        coll="IsoTrackProbe",
                        field=field,
                        bins=21,
                        start=-0.5,
                        stop=20.5,
                        label=f"Probe-track outer {label}",
                    )
                ]
            )
            for field, label in outer_hit_variants.items()
        },
        **{
            f"preMissingOuterTrack_{field}": HistConf(
                [
                    Axis(
                        coll="IsoTrackSearchPreMissingOuter",
                        field=field,
                        bins=21,
                        start=-0.5,
                        stop=20.5,
                        label=f"Search preselection track outer {label}",
                    )
                ]
            )
            for field, label in outer_hit_variants.items()
        },
        **{
            f"preLeptonVetoTrack_{field}": HistConf(
                [
                    Axis(
                        coll="IsoTrackSearchPreLeptonVeto",
                        field=field,
                        bins=21,
                        start=-0.5,
                        stop=20.5,
                        label=f"Search pre-lepton-veto track outer {label}",
                    )
                ]
            )
            for field, label in outer_hit_variants.items()
        },
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
