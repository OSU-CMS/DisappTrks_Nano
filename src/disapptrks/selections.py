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
    """Run-3 TightLepVeto jet ID for NanoAOD jets.

    Prefer the stored NanoAOD ``jetId`` bitmask when available.  Some OSUNano
    productions do not keep the multiplicity branches needed to reconstruct the
    ID from raw fractions, while ``jetId`` is exactly the compact branch meant
    for this use case.  If ``jetId`` is unavailable, fall back to the explicit
    Run-3 TightLepVeto reconstruction used by the first smoke-test files.
    """
    if "jetId" in jets.fields:
        return (jets.jetId & 6) == 6

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


def muon_tag_mask(
    muons, *, pt_min: float = 26.0, eta_max: float = 2.1, iso_max: float = 0.15
):
    """Tight, isolated muon tag used for the first tag-and-probe prototype."""
    return (
        (muons.pt > pt_min)
        & (abs(muons.eta) < eta_max)
        & muons.tightId
        & (muons.pfRelIso04_all < iso_max)
    )


def muon_veto_probe_track_mask(tracks, *, layer: str = "combinedBins"):
    """Probe-track denominator for a first muon-veto tag-and-probe study.

    This intentionally does not apply the muon veto.  The muon-veto pass/fail
    decision is measured on top of this denominator.
    """
    return base_probe_track_mask(
        tracks,
        pt_min=30.0,
        layer=layer,
        apply_calo_cut=True,
        apply_outer_hits_cut=False,
    ) & (
        ((tracks.dRMinElectron < 0.0) | (tracks.dRMinElectron > 0.15))
        & ((tracks.dRMinTauHad < 0.0) | (tracks.dRMinTauHad > 0.15))
    )


def invariant_mass(
    first, second, *, first_mass: float = 0.105658, second_mass: float = 0.105658
):
    """Compute invariant mass from NanoAOD-style pt/eta/phi records."""
    first_px = first.pt * np.cos(first.phi)
    first_py = first.pt * np.sin(first.phi)
    first_pz = first.pt * np.sinh(first.eta)
    first_e = np.sqrt(first_px**2 + first_py**2 + first_pz**2 + first_mass**2)

    second_px = second.pt * np.cos(second.phi)
    second_py = second.pt * np.sin(second.phi)
    second_pz = second.pt * np.sinh(second.eta)
    second_e = np.sqrt(second_px**2 + second_py**2 + second_pz**2 + second_mass**2)

    mass2 = (
        (first_e + second_e) ** 2
        - (first_px + second_px) ** 2
        - (first_py + second_py) ** 2
        - (first_pz + second_pz) ** 2
    )
    return np.sqrt(np.maximum(mass2, 0.0))


def build_muon_veto_tag_probe_pairs(tags, probes):
    """Build flat per-event muon-tag/probe-track pairs for P(muon veto)."""
    import awkward as ak

    tag, probe = ak.unzip(ak.cartesian([tags, probes], axis=1))
    mass = invariant_mass(tag, probe)
    return ak.zip(
        {
            "mass": mass,
            "os": tag.charge * probe.charge < 0,
            "ss": tag.charge * probe.charge > 0,
            "probe_pt": probe.pt,
            "probe_eta": probe.eta,
            "probe_phi": probe.phi,
            "probe_dRMinMuon": probe.dRMinMuon,
            "probe_missingOuterHits": probe.missingOuterHits,
            "probe_passMuonVeto": (probe.dRMinMuon < 0.0) | (probe.dRMinMuon > 0.15),
        }
    )


def os_muon_probe_pair_mask(pairs):
    return pairs.os


def ss_muon_probe_pair_mask(pairs):
    return pairs.ss


def mass10_muon_probe_pair_mask(pairs):
    return pairs.os & (pairs.mass > 10.0)


def ss_mass10_muon_probe_pair_mask(pairs):
    return pairs.ss & (pairs.mass > 10.0)


def z_window_muon_probe_pair_mask(
    pairs, *, z_mass: float = 91.1876, window: float = 10.0
):
    return pairs.os & (pairs.mass > z_mass - window) & (pairs.mass < z_mass + window)


def ss_z_window_muon_probe_pair_mask(
    pairs, *, z_mass: float = 91.1876, window: float = 10.0
):
    return pairs.ss & (pairs.mass > z_mass - window) & (pairs.mass < z_mass + window)


def muon_veto_pair_pass_mask(pairs):
    return pairs.probe_passMuonVeto


def muon_veto_pair_fail_mask(pairs):
    return ~pairs.probe_passMuonVeto


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


def search_track_cutflow_masks(tracks, *, layer: str = "combinedBins"):
    """Return cumulative track masks for debugging the search-track selection."""
    masks = {}
    mask = tracks.pt > 55.0
    masks["track_pt55"] = mask

    mask = mask & (abs(tracks.eta) < 2.1)
    masks["track_eta2p1"] = mask

    mask = mask & ~tracks.inECALCrack
    masks["track_noECALCrack"] = mask

    mask = mask & ~tracks.inDTWheelGap
    masks["track_noDTWheelGap"] = mask

    mask = mask & ~tracks.inCSCTransition
    masks["track_noCSCTransition"] = mask

    mask = mask & ~tracks.inTOBCrack
    masks["track_noTOBCrack"] = mask

    mask = mask & tracks.isFiducialECALTrack
    masks["track_fiducialECAL"] = mask

    mask = mask & (tracks.hp_nValidPixelHits >= 4)
    masks["track_pixelHits4"] = mask

    mask = mask & (tracks.hp_nValidHits >= 4)
    masks["track_validHits4"] = mask

    mask = mask & (tracks.missingInnerHits == 0)
    masks["track_noMissingInner"] = mask

    mask = mask & (tracks.missingMiddleHits == 0)
    masks["track_noMissingMiddle"] = mask

    mask = mask & (tracks.pfRelIso03_chg < 0.05)
    masks["track_chargedIso0p05"] = mask

    mask = mask & (abs(tracks.dxy) < 0.02)
    masks["track_dxy0p02"] = mask

    mask = mask & (abs(tracks.dz) < 0.5)
    masks["track_dz0p5"] = mask

    mask = mask & ((tracks.dRMinJet < 0.0) | (tracks.dRMinJet > 0.5))
    masks["track_dRJet0p5"] = mask

    mask = mask & layer_mask(tracks, layer)
    masks["track_layers4plus"] = mask

    mask = mask & (tracks.caloEnergy < 10.0)
    masks["track_calo10"] = mask

    mask = mask & (tracks.missingOuterHits >= 3)
    masks["track_missingOuter3"] = mask

    mask = mask & ((tracks.dRMinElectron < 0.0) | (tracks.dRMinElectron > 0.15))
    masks["track_electronVeto"] = mask

    mask = mask & ((tracks.dRMinMuon < 0.0) | (tracks.dRMinMuon > 0.15))
    masks["track_muonVeto"] = mask

    mask = mask & ((tracks.dRMinTauHad < 0.0) | (tracks.dRMinTauHad > 0.15))
    masks["track_tauVeto"] = mask

    return masks


def search_event_cutflow_masks(
    analysis_event,
    *,
    met_min: float = 120.0,
    jet_pt_min: float = 110.0,
    jet_met_dphi_min: float = 0.5,
    dijet_dphi_max: float = 2.5,
):
    """Return cumulative event masks for debugging the search event selection."""
    masks = {}

    mask = analysis_event.METNoMu_pt >= met_min
    masks["event_metNoMu120"] = mask

    mask = mask & (analysis_event.leadingJet_pt > jet_pt_min)
    masks["event_leadingJet110"] = mask

    mask = mask & (analysis_event.leadingJetMETNoMuDeltaPhi >= jet_met_dphi_min)
    masks["event_jetMetDphi0p5"] = mask

    mask = mask & (
        (analysis_event.dijetMaxDeltaPhi < 0.0)
        | (analysis_event.dijetMaxDeltaPhi < dijet_dphi_max)
    )
    masks["event_dijetDphi2p5"] = mask

    return masks


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
