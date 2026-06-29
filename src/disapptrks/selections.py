"""Awkward-array selections shared by PocketCoffea and validation tools."""

from __future__ import annotations

import numpy as np


def delta_phi(phi1, phi2):
    return np.arctan2(np.sin(phi1 - phi2), np.cos(phi1 - phi2))


def minimum_delta_r(tracks, objects, object_mask=None):
    """Return the minimum ΔR for every track, using -1 when no object exists."""
    import awkward as ak

    if object_mask is not None:
        objects = objects[object_mask]
    track_pairs, object_pairs = ak.unzip(
        ak.cartesian([tracks, objects], nested=True)
    )
    dr = np.sqrt(
        (track_pairs.eta - object_pairs.eta) ** 2
        + delta_phi(track_pairs.phi, object_pairs.phi) ** 2
    )
    return ak.fill_none(ak.min(dr, axis=2, mask_identity=True), -1.0)


def run3_tight_lepton_veto_jet_mask(jets):
    """Run-3 TightLepVeto jet ID reconstructed from NanoAOD fractions."""
    abs_eta = abs(jets.eta)
    multiplicity = jets.chMultiplicity + jets.neMultiplicity
    central = (
        (abs_eta <= 2.6)
        & (jets.neHEF < 0.99)
        & (jets.neEmEF < 0.90)
        & (multiplicity > 1)
        & (jets.muEF < 0.8)
        & (jets.chHEF > 0.01)
        & (jets.chMultiplicity > 0)
        & (jets.chEmEF < 0.80)
    )
    transition = (
        (abs_eta > 2.6)
        & (abs_eta <= 2.7)
        & (jets.neHEF < 0.90)
        & (jets.neEmEF < 0.99)
        & (jets.muEF < 0.8)
        & (jets.chMultiplicity > 0)
        & (jets.chEmEF < 0.80)
    )
    forward = (
        (abs_eta > 2.7)
        & (abs_eta <= 3.0)
        & (jets.neHEF < 0.99)
        & (jets.neEmEF < 0.99)
        & (jets.neMultiplicity > 1)
    )
    very_forward = (
        (abs_eta > 3.0)
        & (abs_eta <= 5.0)
        & (jets.neEmEF < 0.4)
        & (jets.neMultiplicity > 10)
    )
    return central | transition | forward | very_forward


def layer_mask(tracks, layer: str):
    layers = tracks.hp_nValidTrackerHits
    # Prefer the explicit layer count when supplied by the custom extension.
    if "hp_trackerLayersWithMeasurement" in tracks.fields:
        layers = tracks.hp_trackerLayersWithMeasurement
    if layer == "NLayers4":
        return layers == 4
    if layer == "NLayers5":
        return layers == 5
    if layer == "NLayers6plus":
        return layers >= 6
    if layer == "combinedBins":
        return layers >= 4
    raise ValueError(f"unknown layer bin: {layer}")


def add_isotrack_derived_fields(events):
    """Attach transparent analysis quantities to the ``IsoTrack`` collection."""
    import awkward as ak

    tracks = events.IsoTrack
    abs_eta = abs(tracks.eta)
    tracks = ak.with_field(tracks, (abs_eta >= 1.42) & (abs_eta <= 1.65), "inECALCrack")
    tracks = ak.with_field(tracks, (abs_eta >= 0.15) & (abs_eta <= 0.35), "inDTWheelGap")
    tracks = ak.with_field(tracks, (abs_eta >= 1.55) & (abs_eta <= 1.85), "inCSCTransition")
    tracks = ak.with_field(
        tracks,
        (abs(tracks.dz) < 0.5) & (abs(np.pi / 2.0 - tracks.theta) < 1.0e-3),
        "inTOBCrack",
    )
    tracks = ak.with_field(tracks, tracks.caloEm + tracks.caloHad, "caloEnergy")

    good_jets = (
        (events.Jet.pt > 30.0)
        & (abs(events.Jet.eta) < 4.5)
        & run3_tight_lepton_veto_jet_mask(events.Jet)
    )
    tracks = ak.with_field(
        tracks, minimum_delta_r(tracks, events.Jet, good_jets), "dRMinJet"
    )
    tracks = ak.with_field(
        tracks, minimum_delta_r(tracks, events.Electron), "dRMinElectron"
    )
    tracks = ak.with_field(
        tracks, minimum_delta_r(tracks, events.Muon), "dRMinMuon"
    )
    tracks = ak.with_field(
        tracks, minimum_delta_r(tracks, events.Tau), "dRMinTauHad"
    )
    return tracks


def base_probe_track_mask(
    tracks,
    *,
    pt_min: float = 30.0,
    layer: str = "combinedBins",
    apply_calo_cut: bool = True,
    apply_outer_hits_cut: bool = False,
):
    mask = (
        (tracks.pt > pt_min)
        & (abs(tracks.eta) < 2.1)
        & ~tracks.inECALCrack
        & ~tracks.inDTWheelGap
        & ~tracks.inCSCTransition
        & ~tracks.inTOBCrack
        & tracks.isFiducialECALTrack
        & (tracks.hp_nValidPixelHits >= 4)
        & (tracks.hp_nValidHits >= 4)
        & (tracks.missingInnerHits == 0)
        & (tracks.missingMiddleHits == 0)
        & (tracks.pfRelIso03_chg < 0.05)
        & (abs(tracks.dxy) < 0.02)
        & (abs(tracks.dz) < 0.5)
        & ((tracks.dRMinJet < 0.0) | (tracks.dRMinJet > 0.5))
        & layer_mask(tracks, layer)
    )
    if apply_calo_cut:
        mask = mask & (tracks.caloEnergy < 10.0)
    if apply_outer_hits_cut:
        mask = mask & (tracks.missingOuterHits >= 3)
    return mask


def search_track_mask(tracks, *, layer: str = "combinedBins"):
    return base_probe_track_mask(
        tracks,
        pt_min=55.0,
        layer=layer,
        apply_calo_cut=True,
        apply_outer_hits_cut=True,
    ) & (
        ((tracks.dRMinElectron < 0.0) | (tracks.dRMinElectron > 0.15))
        & ((tracks.dRMinMuon < 0.0) | (tracks.dRMinMuon > 0.15))
        & ((tracks.dRMinTauHad < 0.0) | (tracks.dRMinTauHad > 0.15))
    )


def add_event_derived_fields(events):
    """Build no-muon-MET/jet angular quantities without a custom event table."""
    import awkward as ak

    good = (
        (events.Jet.pt > 30.0)
        & (abs(events.Jet.eta) < 4.5)
        & run3_tight_lepton_veto_jet_mask(events.Jet)
    )
    jets = events.Jet[good]
    order = ak.argsort(jets.pt, ascending=False)
    jets = jets[order]
    leading_pt = ak.fill_none(ak.firsts(jets.pt), -1.0)
    leading_phi = ak.fill_none(ak.firsts(jets.phi), 0.0)

    first, second = ak.unzip(ak.combinations(jets.phi, 2, axis=1))
    dijet_max = ak.fill_none(
        ak.max(abs(delta_phi(first, second)), axis=1, mask_identity=True), -1.0
    )
    return ak.zip(
        {
            "METNoMu_pt": events.MetNoMu.pt,
            "METNoMu_phi": events.MetNoMu.phi,
            "leadingJet_pt": leading_pt,
            "leadingJet_phi": leading_phi,
            "dijetMaxDeltaPhi": dijet_max,
            "leadingJetMETNoMuDeltaPhi": abs(
                delta_phi(leading_phi, events.MetNoMu.phi)
            ),
        }
    )
