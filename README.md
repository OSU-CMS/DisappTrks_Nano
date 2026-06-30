# DisappTrks_Nano

This repository is the migration layer between the Run 3 disappearing-track
analysis and custom NanoAOD/PocketCoffea. The disappearing-track candidates are
the custom NanoAOD `IsoTrack` collection. Object matching, event-level angular
variables, and analysis selections are computed downstream in Python.

## Current status

- The `IsoTrack` branch contract is encoded in `disapptrks.schema`.
- The nominal search-track and event selections are implemented in
  `disapptrks.selections`.
- The legacy fiducial-map hot-spot algorithm is implemented in
  `disapptrks.fiducial`.
- `pocket_coffea` contains the first PocketCoffea workflow and local
  dataset configuration.

The current local smoke-test file is expected outside this repository as
`../nano_99.root`. Its dataset metadata is labeled as 2024F Muon data. Check
productions before running:

```bash
python3 -m pip install -e '.[analysis,test,pocket-coffea]'
disapptrks audit-schema ../nano_99.root
```

The audit exits with status 2 if the requested scope is incomplete.

Build a PocketCoffea dataset JSON directly from an EOS directory:

```bash
disapptrks make-dataset-json /store/user/YOUR_PATH \
  --recursive \
  --dataset-name Run2024G_Muon_OSUNano_EOS \
  --sample DATA_Muon \
  --year 2024 \
  --era G \
  --primary-dataset Muon \
  -o pocket_coffea/datasets/eos_2024G_muon.json
```

The command uses `xrdfs root://cmseos.fnal.gov ls -u` so the output JSON gets
full XRootD file URLs. If you already have a text filelist, pass it with
`--filelist` instead of an EOS path.

## PocketCoffea setup

Install this analysis and the recommended PocketCoffea release in one
environment:

```bash
python3 -m pip install -e '.[analysis,test,pocket-coffea]'
```

The `pocket-coffea` extra pins `pocket-coffea==0.9.12`, matching the current
PocketCoffea release line and its Coffea/Awkward dependency stack. If you are
developing PocketCoffea itself from a sibling checkout, install that checkout
instead of the extra:

```bash
python3 -m pip install -e ../PocketCoffea
python3 -m pip install -e '.[analysis,test]'
```

Once a file passes the schema audit:

```bash
cd pocket_coffea
pocket-coffea run --cfg config.py --test -lf 1 -lc 1 -c 50000 -e iterative -ps -o output_test_nano99
```

Summarize the current muon-veto tag-and-probe prototype from a PocketCoffea
output file:

```bash
disapptrks summarize-pveto output_test_nano99/output_Run2024F_Muon_OSUNano_local.coffea
```

By default this reports
`muon_veto_zwindow_pass / muon_veto_zwindow` with a binomial statistical
uncertainty.

To apply the same-sign subtraction used to estimate the flat non-Z background
under the opposite-sign Z peak:

```bash
disapptrks summarize-pveto output_test_nano99/output_Run2024F_Muon_OSUNano_local.coffea --ss-subtract
```

This reports
`(muon_veto_zwindow_pass - muon_veto_ss_zwindow_pass) / (muon_veto_zwindow - muon_veto_ss_zwindow)`.

## Physics-equivalence caveat

For data, the observed `IsoTrack` missing-hit values can be used directly.
Exact MC equivalence still requires deciding how to reproduce the stochastic
hit-drop and TOB-drop treatment from the legacy OSUT3 workflow. The electron
and muon inefficiency maps are analysis inputs applied downstream; they are not
assumed to be embedded in a separate track table.

See [MIGRATION.md](MIGRATION.md) for the source-to-PocketCoffea mapping and
validation sequence.

Run unit tests with:

```bash
python3 -m pip install -e '.[test]'
pytest
```
