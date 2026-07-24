# Codex Handoff: Background Estimate Development

This note is intended to be handed to another Codex session before continuing
development of the Run-3 DisappTrks PocketCoffea background-estimate workflow.
It summarizes the code structure, conventions, current modes, and known
debugging points.

Companion files:

- `docs/codex_memory_background_estimates.md`: compact project memory/state.
- `docs/codex_skills_background_estimates.md`: operating guidance for a Codex
  session working on this code.

## Repository Context

There are two related codebases under the same parent directory:

- `DisappTrks/`: the legacy CMSSW/OSUT3 DisappTrks analysis code.
- `DisappTrks_Nano/`: the newer custom NanoAOD plus PocketCoffea migration.

When implementing or debugging the Nano workflow, use the legacy code as the
physics reference, especially:

- `DisappTrks/BackgroundEstimation/python/`
- `DisappTrks/BackgroundEstimation/test/bkgdEstimate_2022.py`
- `DisappTrks/StandardAnalysis/python/EventSelections.py`
- `DisappTrks/StandardAnalysis/python/Cuts.py`

The Nano code should keep AN terminology in public names where possible:
`Pveto`, `Poffline`, `Pmiss`, `basic_selection`,
`isolated_track_selection`, `candidate_track_selection`, and
`disappearing_track_selection`.

## Main Nano Files

| File | Role |
| --- | --- |
| `src/disapptrks/selections.py` | Physics definitions on Awkward arrays: object masks, tag definitions, probe-track masks, pair builders, AN-style cutflow masks. |
| `pocket_coffea/workflow.py` | Builds derived event/object collections and event-level counters during PocketCoffea processing. |
| `pocket_coffea/cuts.py` | PocketCoffea `Cut` objects, usually thin wrappers around counts or booleans stored by `workflow.py`. |
| `pocket_coffea/config.py` | Dataset filtering, category mode routing, category selection, histogram variable selection. |
| `src/disapptrks/lepton_backgrounds.py` | Computes final lepton background estimates from `Pveto`, `Poffline`, and `Pmiss` inputs. |
| `src/disapptrks/tables.py` | Cutflow and Pveto table extraction/formatting. |
| `src/disapptrks/cli.py` | CLI commands for fiducial maps, Pveto tables, fake-track tables, and lepton-background estimates. |
| `docs/pocket_coffea_workflows.md` | Collaborator-facing workflow guide. Keep it updated when modes or commands change. |

Rule of thumb:

- Put reusable physics masks in `selections.py`.
- Put derived collections and counters in `workflow.py`.
- Put named PocketCoffea categories in `cuts.py` and `config.py`.
- Put postprocessing formulas in `lepton_backgrounds.py` or `tables.py`.

## Current Production Modes

The PocketCoffea mode is selected with:

```bash
DISAPPTRKS_CATEGORY_MODE=<mode>
```

Important modes:

| Mode | Dataset | Purpose |
| --- | --- | --- |
| `fiducial_maps` | `DATA_Muon` or `DATA_EGamma` | Produce before/after eta-phi histograms used to build electron and muon fiducial-map JSON/NPZ files. |
| `muon_pveto` | `DATA_Muon` | Muon `Pveto` tag-probe pairs and categories. |
| `electron_pveto` | `DATA_EGamma` | Electron `Pveto` tag-probe pairs and categories. |
| `tau_mu_pveto` | `DATA_Muon` | Tau `Pveto` measurement using muon low-`MT` tags. |
| `tau_ele_pveto` | `DATA_EGamma` | Tau `Pveto` measurement using electron low-`MT` tags. |
| `muon_pmiss_poffline` | `DATA_Muon` | Muon `Poffline` and `Pmiss` control categories only. |
| `electron_pmiss_poffline` | `DATA_EGamma` | Electron `Poffline` and `Pmiss` control categories only. |
| `tau_mu_pmiss_poffline` | `DATA_Muon` | Tau control `Poffline` and `Pmiss` with muon low-`MT` tags. |
| `tau_ele_pmiss_poffline` | `DATA_EGamma` | Tau control `Poffline` and `Pmiss` with electron low-`MT` tags. |
| `fake_tracks` | `DATA_JetMET`, `DATA_MET`, `DATA_Muon`, or `DATA_EGamma` | Fake-track estimate control regions. |

Avoid using `muon_backgrounds`, `egamma_backgrounds`, or `all` for production
unless explicitly debugging. They can overload PocketCoffea/Coffea because they
select too many categories and can exhaust `PackedSelection` slots.

## Fiducial Maps

The fiducial-map job stores:

- `ElectronFiducialBefore`
- `ElectronFiducialAfter`
- `MuonFiducialBefore`
- `MuonFiducialAfter`

The JSON/NPZ maps are built from the `.coffea` output with:

```bash
disapptrks make-fiducial-map \
  --flavor electron \
  --output-json /path/to/electron_fiducial_map.json \
  /path/to/fiducial_map_output/output_*.coffea
```

or `--flavor muon`.

Hot spots use the default threshold:

```bash
--threshold 2.0
```

For eras with one pathological high-inefficiency bin inflating the stddev, use:

```bash
--stddev-exclude-top 1
```

This removes the most extreme occupied bin only from the stddev calculation.
The bin is still tested and reported as a hot spot.

Pveto jobs load maps with either explicit paths:

```bash
DISAPPTRKS_ELECTRON_FIDUCIAL_MAP_JSON=/path/to/electron_fiducial_map.json
DISAPPTRKS_MUON_FIDUCIAL_MAP_JSON=/path/to/muon_fiducial_map.json
```

or a directory:

```bash
DISAPPTRKS_FIDUCIAL_MAP_DIR=/path/to/fiducial_maps
```

The directory must contain exactly:

```text
electron_fiducial_map.json
muon_fiducial_map.json
```

For production validation, set:

```bash
DISAPPTRKS_REQUIRE_FIDUCIAL_MAPS=1
```

This makes the job fail if the maps cannot be loaded or contain no hot spots.

## Lepton Background Formula

The AN-style estimate is:

```text
N_lepton = N_ctrl * Pveto * Poffline * Pmiss / epsilon_trig^lepton
```

In Nano postprocessing this is handled by:

```text
src/disapptrks/lepton_backgrounds.py
```

Definitions:

- `N_ctrl`: control yield category, e.g.
  `electron_background_control_NLayers4`.
- `Poffline`: ratio of offline-MET-passing control events to control events.
- `Pmiss`: ratio of MET-trigger-passing control events to offline-MET-passing
  control events.
- `Pveto`: probability from the tag-probe pair categories.
- `epsilon_trig^lepton`: lepton trigger efficiency. The postprocessor default
  is `1.0`; pass the AN/legacy value with `--trigger-efficiency` and
  `--trigger-efficiency-error`.

`N_ctrl` can also be scaled with `--control-prescale`. This matches the legacy
MET/lepton-dataset luminosity or prescale correction. The default is `1.0`.

Current category naming:

```text
<mode>_background_control_{layer}
<mode>_background_offline_{layer}
<mode>_background_trigger_{layer}
```

where `<mode>` is `muon`, `electron`, `tau_mu`, or `tau_ele`, and `{layer}` is
`NLayers4`, `NLayers5`, `NLayers6plus`, or `combinedBins`.

## Important Legacy Convention For Pveto

For electron and muon backgrounds, the legacy code does not use the direct
tag-probe ratio for `Pveto`. It uses:

```text
Pveto = N_pass / (2*N_total - N_pass)
```

after same-sign subtraction.

Reference:

```text
DisappTrks/BackgroundEstimation/python/bkgdEstimate.py
LeptonBkgdEstimate.printPpassVetoTagProbe()
```

The relevant legacy branch is:

```python
if (self._flavor == "electron" or self._flavor == "muon") and not self._useHistogramsForPpassVeto:
    eff = scaledPasses / (2.0 * total - scaledPasses)
else:
    eff = scaledPasses / total
```

Tau uses the direct `N_pass / N_total` branch.

Nano implements this in `src/disapptrks/lepton_backgrounds.py` via
`pveto_count_from_pair_counts(..., use_two_lepton_denominator=True)` for
electron/muon flavors.

## Poffline/Pmiss Control Selection

The control-only modes build the relevant tag collection and count events with
near-lepton tracks. The core helper is:

```text
pocket_coffea/workflow.py::_lepton_background_track_mask
```

Current behavior:

- track `pt > 55 GeV`
- requested layer bin
- standard isolated-track quality
- `dR(track, jet) > 0.5`
- no missing-outer-hit requirement
- for muons: `caloEnergy < 10 GeV`
- for electrons: no calorimeter-energy cut
- electron control: `0 <= dRMinElectron < 0.1`
- muon control: `0 <= dRMinMuon < 0.1`
- tau control: `0 <= dRMinTauHad < 0.1`

The `dR(track, jet) > 0.5` requirement is important. It was added because the
legacy `ElectronTagPt55` and `MuonTagPt55` control channels include
`isoTrkCuts`, which include `cutTrkJetDeltaPhi`.

`Poffline` and `Pmiss` are counted in:

```text
pocket_coffea/workflow.py::_store_lepton_background_controls
```

The offline event uses:

- event quality: golden JSON, MET filters, jet-veto map
- at least one tag
- MET-no-mu-minus-selected-lepton `pt >= 120 GeV`
- leading-jet delta phi to that MET direction `>= 0.5`

The trigger event is:

```text
offline_event AND _met_trigger_mask(events)
```

New outputs from the `*_pmiss_poffline` modes also store the histograms needed
to reproduce the legacy integration:

```text
n<Prefix>BackgroundMetMinusOnePt_{layer}
n<Prefix>BackgroundMetMinusOnePtTrig_{layer}
n<Prefix>BackgroundDeltaPhiMetJetLeadingVsMetMinusOnePt_{layer}
```

where `<Prefix>` is `Muon`, `Electron`, `TauMu`, or `TauEle`.

The postprocessor uses these histograms automatically when present. It builds a
MET-trigger turn-on from the first two histograms, weights the 2D
lepton-removed MET versus delta-phi histogram, and integrates the region above
`--met-cut` and `--phi-cut`. If the histograms are missing, it falls back to
the older scalar cutflow ratios and prints `met_method=cutflow-ratio`.

After changing this code, rerun the relevant `*_pmiss_poffline` jobs. Existing
Pveto outputs can be reused.

For 2022 electrons, the legacy script calls `useFilesForTriggerEfficiency()`
unless the flat trigger-efficiency option is enabled. The flat fallback in that
script is `0.840 +/- 0.005`, but the file-derived value is the nominal AN-style
choice.

## Typical Commands

Electron Pveto:

```bash
cd DisappTrks_Nano/pocket_coffea
DISAPPTRKS_CATEGORY_MODE=electron_pveto \
DISAPPTRKS_ENABLE_PVETO_DIAGNOSTICS=0 \
DISAPPTRKS_REQUIRE_FIDUCIAL_MAPS=1 \
DISAPPTRKS_DATASET_JSON=datasets/eos_2022CD_EGamma.json \
DISAPPTRKS_FIDUCIAL_MAP_DIR=/path/to/fiducial_maps/2022CD \
python -m pocket_coffea.scripts.runner run \
  --cfg config.py \
  --outputdir analysis_output/2022CD_electron_pveto \
  --executor dask@lpc \
  --executor-custom-setup executors_lpc.py \
  --custom-run-options run_options_lpc_dask.yaml \
  --scaleout 60 \
  --queue workday
```

Electron Poffline/Pmiss:

```bash
cd DisappTrks_Nano/pocket_coffea
DISAPPTRKS_CATEGORY_MODE=electron_pmiss_poffline \
DISAPPTRKS_ENABLE_PVETO_DIAGNOSTICS=0 \
DISAPPTRKS_DATASET_JSON=datasets/eos_2022CD_EGamma.json \
python -m pocket_coffea.scripts.runner run \
  --cfg config.py \
  --outputdir analysis_output/2022CD_electron_pmiss_poffline \
  --executor dask@lpc \
  --executor-custom-setup executors_lpc.py \
  --custom-run-options run_options_lpc_dask.yaml \
  --scaleout 60 \
  --queue workday
```

Postprocess electron background:

```bash
cd DisappTrks_Nano
disapptrks estimate-lepton-background \
  --mode electron \
  --run-period 2022CD \
  --trigger-efficiency <epsilon> \
  --trigger-efficiency-error <epsilon_error> \
  --output-json tables/electron_background_2022CD.json \
  --output-tex tables/electron_background_2022CD.tex \
  pocket_coffea/analysis_output/2022CD_electron_pveto/output_*.coffea \
  pocket_coffea/analysis_output/2022CD_electron_pmiss_poffline/output_*.coffea
```

For a smoke test, add this to the runner command:

```bash
--limit-files 1 --limit-chunks 1 --scaleout 2 --queue microcentury
```

## Debugging Checklist

If a job fails with:

```text
RuntimeError: Exhausted all slots in PackedSelection
```

then the selected category set is too large. Use a focused mode such as
`electron_pveto` or `electron_pmiss_poffline`, and keep diagnostics disabled:

```bash
DISAPPTRKS_ENABLE_PVETO_DIAGNOSTICS=0
DISAPPTRKS_ENABLE_SEARCH_DIAGNOSTICS=0
```

If histogram filling fails with a missing key like:

```text
ValueError: key "nElectronFiducialHotSpotsLoaded" does not exist
```

then `config.py` selected a histogram variable that `workflow.py` did not
create in that mode. Fix either the variable filtering in `_variables_for_mode`
or create the event field in `apply_object_preselection`/`count_objects`.

If the fiducial map seems to have no effect:

- Run with `DISAPPTRKS_REQUIRE_FIDUCIAL_MAPS=1`.
- Check `nElectronFiducialHotSpotsLoaded` and
  `nMuonFiducialHotSpotsLoaded`.
- For muon Pveto, compare:
  `nMuonPVetoTagProbePairZWindowPassNoFiducial` and
  `nMuonPVetoTagProbePairZWindowFiducialRejected`.

If electron 2022CD numbers do not match AN Table 28:

- Check the `Pveto` convention first:
  electron/muon should use `N_pass/(2*N_total - N_pass)` after SS subtraction.
- Check that the Poffline/Pmiss control track includes `dR(track, jet) > 0.5`.
- Check the legacy prescale/luminosity factor:
  `MET lumi / EGamma lumi` is used in legacy 2022 electron estimates.
- Check whether legacy `Pmiss` used trigger-efficiency files rather than a
  direct MET-HLT event bit.
- Compare the exact electron tag definition:
  Nano uses `Electron.cutBased >= 4`, dxy/dz barrel/endcap cuts, `pt > 35`,
  `|eta| < 2.1`, and the single-electron HLT mask.

## Development Safety

The working tree may contain unrelated local changes. Before editing, run:

```bash
git -C DisappTrks_Nano status --short
```

Do not revert unrelated files. Keep edits tightly scoped, and use:

```bash
python -m py_compile <changed .py files>
python -m pytest DisappTrks_Nano/tests/test_lepton_backgrounds.py
git -C DisappTrks_Nano diff --check
```

for quick validation.

## Current Open Validation Items

The Nano implementation now has the basic plumbing for fiducial maps, `Pveto`,
`Poffline`, and `Pmiss`, but the following should still be validated against
the AN and/or legacy outputs:

- Electron 2022CD Table-28 closure after rerunning
  `electron_pmiss_poffline`.
- Whether the Nano `Pmiss` direct trigger-bit ratio needs to be replaced by the
  legacy trigger-efficiency-file method.
- Whether the legacy `MET lumi / lepton lumi` prescale factor should be applied
  explicitly in Nano postprocessing.
- The exact Run-3 electron veto object used for fiducial maps:
  Nano currently uses `Electron.cutBased >= 1` for veto electrons.
- The exact Run-3 tau object working point in Nano:
  `hadronic_tau_veto_object_mask` uses DeepTau 2018v2p5 raw thresholds and
  `idDecayModeNewDMs`.
