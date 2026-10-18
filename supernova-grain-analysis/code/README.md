# Supernova particle data audit

This project maps the requested 50AM, jet3b, cco2 and vconvL data to their actual
files on ASU SOL. “Grain” means an **SPH simulation particle**, as clarified by
the project owner. The audit preserves the observed `vconv` label until its
association with vconvL is established for the relevant product.

**Start with the local visual report:**
`/scratch/rpraka15/supernova-grain-analysis/reports/index.html`.
Open it in a browser on SOL, or download the project scratch `reports`, `plots`
and `tables` directories together. It is self-contained apart from those local
files and does not contact external services.

## What was found

- The actual explosions have roughly 0.93–0.95 million particles per selected
  snapshot. The old `50AM_test` reports inspected a bundled Noh test instead.
- SDF snapshots supply the hydrodynamic fields. A numeric PRISM directory
  supplies one particle's postprocessed composition, with `x.txt` for isotope
  mass fractions and `xz.txt` for elemental mass fractions.
- `xtime.txt` supplies time-resolved composition and saved thermodynamic states.
  Final outputs are at different configured times, not a common physical epoch.
- Your 298 requested fields are mapped in
  `tables/requested_field_map.csv` in scratch. The sampled values table retains
  missingness and source paths. It is **not an all-particle production dataset**.
- Current trajectory scripts and sampled inputs show serious time/density
  discrepancies. These could affect postprocessed abundances; fixing a label
  after nucleosynthesis does not repair a calculation run with wrong inputs.
- Stored late-time `Y_el` and hydro composition normalization also need checks.
  Never silently replace a stored zero, renormalize, or splice different epochs.

## Read the evidence

- [Hydro provenance, layouts and unit checks](docs/provenance.md)
- [PRISM formats, sparse rows and older CSV pitfalls](docs/prism-formats.md)
- [Particle census and archive checks](docs/coverage.md)
- [Audit interpretation and extraction contract](docs/analysis-guide.md)

The audit is deliberately read-only with respect to source data. No collaborator
scripts are executed. No files are deleted. Code and documentation live here;
generated data, tables, logs, plots and a code snapshot live in your scratch
project. Raw simulations and archives are never copied into Git.

## Reproduce on SOL

The existing environment is `/home/rpraka15/.conda/envs/ccsn/bin/python`.
These commands run bounded inspection and diagnostic generation, not simulations.

```bash
cd /home/rpraka15/supernova-grain-analysis
/home/rpraka15/.conda/envs/ccsn/bin/python -m unittest discover -s tests -v
/home/rpraka15/.conda/envs/ccsn/bin/python scripts/inventory.py
/home/rpraka15/.conda/envs/ccsn/bin/python scripts/audit.py
/home/rpraka15/.conda/envs/ccsn/bin/python scripts/unit_crosscheck.py
/home/rpraka15/.conda/envs/ccsn/bin/python scripts/build_report.py
```

`audit.py` reads every header in the named snapshot directories, validates
record lengths, checks every ID in one selected snapshot per model, then reads
bounded PRISM samples. Archive reads stop at a 64 MiB decompressed prefix budget
per selected archive and never extract paths. The inventory script separately
lists direct particle entry names; directory counts do not prove job completion.
See each script's `--help` where available. Outputs may be regenerated in place;
source inputs remain untouched.

Runtime versions are recorded in `debug/run_summary.json`; dependencies are in
[requirements.txt](requirements.txt). Tests exercise binary endianness, packed
records, corruption detection and sparse/time-block parsing.

## GitHub status

This is a local Git repository on branch `main`. GitHub publication was deferred
by the owner until authentication. No remote or authenticated account has been
provided; a GitHub push cannot succeed yet. The GitHub CLI was not on PATH during
inspection. No credentials are stored here.

Git author identity was absent. Local commits use repository-only
`rpraka15 <rpraka15@localhost>` until an appropriate publication identity is
provided. This does not change global Git configuration.

When the owner is ready, authenticate through GitHub's browser/device flow,
choose the destination and visibility, set the publication identity and remote,
then push `main`. Do not paste a token into this repository or a chat message.
