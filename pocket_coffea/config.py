"""Minimal PocketCoffea configuration for local custom-NanoAOD validation.

Run from this directory after installing this package and the sibling
PocketCoffea checkout.
"""

from __future__ import annotations

import os

import cloudpickle

from pocket_coffea.parameters import defaults
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.utils.configurator import Configurator

import cuts
import workflow
from cuts import has_disappearing_track, search_kinematics
from workflow import DisappTrksProcessor

cloudpickle.register_pickle_by_value(cuts)
cloudpickle.register_pickle_by_value(workflow)

localdir = os.path.dirname(os.path.abspath(__file__))
parameters = defaults.get_default_parameters()

cfg = Configurator(
    parameters=parameters,
    datasets={
        "jsons": [f"{localdir}/datasets/local_2022D_muon.json"],
        "filter": {"samples": ["DATA_Muon"], "year": ["2022"]},
    },
    workflow=DisappTrksProcessor,
    calibrators=[],
    skim=[],
    preselections=[],
    categories={
        "inclusive": [passthrough],
        "search": [search_kinematics, has_disappearing_track],
    },
    weights={"common": {"inclusive": []}, "bysample": {}},
    weights_classes=[],
    variations={"weights": {"common": {"inclusive": []}}},
    variables={},
    columns={},
)
