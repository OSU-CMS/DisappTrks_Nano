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

## LPC Dask submission

The preferred LPC path is PocketCoffea's Dask executor with the local LPC
executor from `pocket_coffea/executors_lpc.py`.  This follows the working
`displaced_leptons` setup and avoids maintaining one custom Condor job wrapper
per file slice.

From the repository root on `cmslpc`, run the bootstrap once:

```bash
./setup_lpc.sh
./shell
```

Inside `./shell`, install this analysis if needed and run from `pocket_coffea`.
Use `python -m pocket_coffea.scripts.runner run` rather than relying on the
`pocket-coffea` executable; on the LPC container the executable can point at a
Python outside the container/venv and miss `lpcjobqueue`.

```bash
cd pocket_coffea
DISAPPTRKS_DATASET_JSON=datasets/eos_2023C_muon.json \
python -m pocket_coffea.scripts.runner run \
  --cfg config.py \
  --outputdir analysis_output/2023C_muon_pveto \
  --executor dask@lpc \
  --executor-custom-setup executors_lpc.py \
  --custom-run-options run_options_lpc_dask.yaml \
  --limit-files 1 \
  --limit-chunks 1
```

Tune `scaleout`, `chunksize`, memory, and queue defaults in
`pocket_coffea/run_options_lpc_dask.yaml`.  Condor logs are written under
`$HOME/pocketcoffea_dask_logs/<output-tag>/condor_log` by default so they are
visible to the LPC schedd outside the container.

After the smoke test starts workers successfully, increase worker count from the
command line, for example `--scaleout 60 --queue workday`.

## LPC manual Condor fallback

If Dask/lpcjobqueue is unavailable, the v2-style wrapper in
`pocket_coffea/scripts` remains as a fallback.  The Condor worker runs in a
sandbox, so it does not depend on `/uscms` or `/uscms_data` paths being mounted
inside the job.

Build a relocatable Python target directory once from the repository root:

```bash
python3 -m pip install --target pocket_coffea/python_env '.[analysis,pocket-coffea]'
```

If the target environment already exists and EOS reads fail with
`No module named 'XRootD'`, add the XRootD Python bindings to it:

```bash
python3 -m pip install --target pocket_coffea/python_env --upgrade xrootd
```

Then submit from `pocket_coffea`:

```bash
cd pocket_coffea
FILES_PER_JOB=5 CHUNKSIZE=50000 \
  scripts/submit_lpc_dataset.sh datasets/eos_2023C_muon.json 2023C_muon_pveto DATA_Muon 2023
```

The wrapper transfers `config.py`, `cuts.py`, `workflow.py`, `../src`, the
dataset JSON, your X509 proxy, and `python_env` into the Condor sandbox.  Each
job makes a per-job dataset JSON, runs PocketCoffea iterative over its file
slice, and returns outputs under `analysis_output/<tag>/`.

After the jobs finish, merge the per-job outputs:

```bash
pocket-coffea merge-outputs analysis_output/2023C_muon_pveto/*.coffea \
  -o analysis_output/2023C_muon_pveto/output_Run2023C_Muon_OSUNano_EOS.coffea
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
