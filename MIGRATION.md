# Run 3 migration map

## Authoritative legacy sources

The migration currently treats these files as the executable reference:

- `DisappTrks/BackgroundEstimation/test/bkgdEstimate_2023.py`: composition of
  each estimate.
- `DisappTrks/BackgroundEstimation/python/ElectronTagProbeSelections.py`
- `DisappTrks/BackgroundEstimation/python/MuonTagProbeSelections.py`
- `DisappTrks/BackgroundEstimation/python/TauTagProbeSelections.py`
- `DisappTrks/BackgroundEstimation/python/bkgdEstimate.py`: probability
  combination and same-sign subtraction.
- `DisappTrks_v2/BkgdEstimation/scripts/*Background*_v2*.py`: current
  flat-ntuple implementations and cutflows.
- `DisappTrks_v2/BkgdEstimation/scripts/make_fiducial_maps_v2.py`
- `DisappTrks_v2/BkgdEstimation/scripts/fiducial_map_tools.py`

The analysis note should be used as the final arbiter whenever one of these
implementations disagrees with the documented selection.

## Estimate structure

For each lepton flavor and tracker-layer bin, the legacy estimate is

```text
N_est = N_ctrl × scale × P_veto × P_offlineMET × P_METtrigger
```

with the configured trigger-efficiency divisor and era-specific corrections.

- Electron: `Z → ee` tag-and-probe supplies `P_veto`; EGamma control data
  supplies `N_ctrl` and the offline-MET probability.
- Muon: `Z → μμ` tag-and-probe supplies `P_veto`; Muon control data supplies
  `N_ctrl` and the offline-MET probability.
- Tau: the numerator combines `Z → ττ → μ+track` and
  `Z → ττ → e+track` channels, including same-sign subtraction and the
  channel-dependent scale-factor combination in `printPpassVetoTagProbe`.
  Tau-trigger and muon+tau-trigger corrections are additional Run 3 inputs.

The PocketCoffea implementation must retain pair multiplicity. Counting only
events with at least one tag-probe pair is not equivalent to the legacy
histogram integral.

## Fiducial maps

For electron and muon maps:

1. Fill `(eta, phi)` before and after the corresponding lepton veto.
2. Compute `after / before` in occupied bins.
3. Compute the global mean as `sum(after) / sum(before)`.
4. Compute the sample standard deviation across occupied bins.
5. Mark bins more than 2σ above the mean inefficiency.
6. Veto tracks within at least ΔR 0.05 (or the bin half-diagonal if larger).

The numerical implementation is in `src/disapptrks/fiducial.py`.

## Validation gates

1. A NanoAOD file must pass `disapptrks audit-schema`.
2. Reproduce object multiplicities and every cumulative cut row on the same
   MiniAOD events.
3. Compare tag-probe pair counts, including OS and SS pairs, not only events.
4. Compare `beforeVeto` and `afterVeto` map histograms bin by bin.
5. Compare `N_ctrl`, each probability, `alpha`, and `N_est` independently for
   4-layer, 5-layer, 6+-layer, and combined categories.
6. Only then compare the final summed background.

## Current implementation notes

- Candidate tracks are read from `events.IsoTrack`.
- `MetNoMu_*` is used when stored. The schema also permits reconstructing it
  from a base MET collection and the corresponding muon momenta.
- Track-to-jet/electron/muon/tau distances and jet/MET angular variables are
  computed in the PocketCoffea workflow.
- Electron and muon fiducial maps are downstream analysis inputs.
- Exact MC equivalence still needs the legacy hit-drop/TOB-drop treatment.
- Trigger-object matching definitions still need an era-by-era mapping from
  `TrigObj_filterBits` to the legacy tag requirements.
