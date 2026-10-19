# Supernova analysis reference

This repository collects the code and small validation records for
the supernova particle audit, the 50Am chunk exports, and the independent table
comparison performed on 18 October 2026. It is separate from
`ccsn-presolar-oxide-analysis`. Its GitHub destination is the private repository
[rishantp22/supernova-analysis-reference](https://github.com/rishantp22/supernova-analysis-reference).

## Contents

| Directory | Purpose |
|---|---|
| `supernova-grain-analysis/code/` | SDF/PRISM readers, inventory, data audit, unit checks, report generation |
| `supernova-grain-analysis/reports/` | Interpretation, coverage and provenance documentation |
| `supernova-grain-analysis/debug/`, `inventory/`, `tables/` | Selected small diagnostic records and field maps |
| `final_chunks_export_20261018/` | Export, join validation, plotting conversion and packaging scripts, plus small validation summaries |
| `notebooks/` | Two original analysis notebooks with outputs and execution counts cleared in the copies |
| `comparison/` | Streaming comparator, correctness tests, Slurm script, software versions and measured comparison |
| `slurm/` | Portable comparison-job template and original historical copy-job script |
| `schema/` | The requested 298-column schema |
| [`50AM_full_nucleosynthesis_comparison/`](50AM_full_nucleosynthesis_comparison/) | Completed all-particle, all-isotope audit: code, reports, all figures, statistics and provenance |
| `provenance/` | Source hashes, preparation metadata and the complete excluded-file inventory |

The source scripts and small reports are exact copies. Notebook code cells are
unchanged. Original files and datasets remain at their existing scratch paths;
none was moved, deleted, repaired, or overwritten. Historical reports describe
the state when they were written, including old paths and GitHub status. Current
status and current-file discrepancies are documented here and in
[the comparison report](comparison/REPORT.md).

## Reproduce

Use Python 3.10 or later with the packages in `requirements.txt`. The SOL
environment used for validation was `/home/rpraka15/.conda/envs/ccsn/bin/python`;
exact installed versions are recorded in `comparison/environment.json`.

Run the small comparator checks without loading simulation datasets:

```bash
python -m unittest discover -s comparison -p 'test_compare_outputs.py' -v
PYTHONPATH=supernova-grain-analysis/code python -m unittest discover -s tests -v
```

Run the full comparison in Slurm, from this repository directory:

```bash
sbatch slurm/compare_tables.sbatch \
  /scratch/rpraka15/final_chunks.out \
  /scratch/rpraka15/50Am_plotting_all_isotopes.out \
  /scratch/rpraka15/comparison-rerun
```

The output directory must be new. Set `PYTHON_BIN` to the intended Python
executable before submitting if `python3` is not the configured environment.
Supply the site's account/partition options to `sbatch` as needed. The comparison
reads originals only, hashes all bytes, validates records, sorts temporary copies
by integer particle ID, and compares columns by name in bounded memory.

The historical audit/export scripts retain their original SOL paths and input
assumptions. Their `ROOT`, `WORK`, output locations and the paths documented in
their READMEs must be reviewed before rerunning in another environment. They are
reference calculations, not an automatically portable installed package. In
particular, export scripts create new products and some report scripts regenerate
derived files. Do not run them as an installation step. Use a separate scratch
working directory for a new computation and keep generated products outside Git.

Export sequence: build the last-category table with `build_final_chunks.py`,
verify it with `check_final_export.py` and `verify_join_logic.py`, create the
CGS plotting form with `make_plotting_format.py`, and optionally add complete-model
coverage with `add_complete_model_ids.py`. Packaging scripts are preserved for
reference. The historical Slurm copy script is included as provenance and was
not executed during this preparation.

## Data and publication boundaries

Large exports, chunk datasets, per-chunk intermediates, numerical caches,
full particle catalogs, embedded HTML, plotting caches, archives, and credentials are
excluded. `provenance/excluded_files.csv` lists the inspected excluded files and
their sizes. `provenance/source_manifest.csv` maps each copied artifact to its
original path and source/prepared SHA-256. These records identify data dependencies;
they do not bundle the datasets or assert redistribution rights.

Both inspected scratch directories were outside Git. The earlier local repository
`/home/rpraka15/supernova-grain-analysis` was on `main`, clean, and had no remote.
No excluded research dataset was tracked in those inspected source repositories.
The repository uses `main`. The preparation and Git preflight records in
`provenance/` are historical snapshots taken before repository creation and
publication; their upload and remote fields describe that earlier stage.

`NA` is missing information, never zero. Particle IDs refer to simulation/SPH
particles; they do not establish condensed dust grains. Matching IDs does not
make hydrodynamic and nucleosynthesis epochs equivalent. Field definitions and
native abundance epochs are retained in the source documentation.

## Completed full-population audit

The new [50AM audit](50AM_full_nucleosynthesis_comparison/README.md) compares 937,227 matched valid particle IDs, audits 1,876 chunk files, and records the complete available species inventory. Start with its [key findings](50AM_full_nucleosynthesis_comparison/reports/KEY_FINDINGS.md). All 113 PNG/PDF figure pairs are included. Large databases remain in the separately prepared complete export; the [publication manifest](50AM_full_nucleosynthesis_comparison/PUBLICATION_MANIFEST.csv) and [exclusion list](50AM_full_nucleosynthesis_comparison/EXCLUDED_FILES.csv) document the exact boundary.
