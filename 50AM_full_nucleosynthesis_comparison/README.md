# Reproducible 50AM nucleosynthesis comparison

> GitHub contains code, reports and small tables. See [publication scope](PUBLICATION.md) for the complete export, omitted databases and reproduction requirements.

Read [the scientific report](reports/50AM_nucleosynthesis_comparison_report.md) and [key findings](reports/KEY_FINDINGS.md) before using any abundance or yield. These are audits of existing calculations, not newly calculated nucleosynthesis. Both plotting files have unresolved reliability problems.

Original inputs are opened read-only. Scripts, tables, extracted provenance, reports and PNG/PDF figures remain in this directory. No original dataset is repaired, normalized, moved or deleted.

| Directory | Contents |
| --- | --- |
| `scripts/` | Python programs, meaningful parser/numerical tests, Slurm submission files |
| `tables/` | Complete ID-aligned arrays, all-species statistics, yields, ratios, source validation |
| `figures/` | Population density plots and diagnostic figures, PNG and vector PDF |
| `provenance/` | Input fingerprints, complete group metadata inventory, source snapshots, run controls |
| `reports/` | Scientific, provenance, chunk integrity, numerical validation and interpretation reports |
| `logs/` | Slurm stdout/stderr, including preserved unsuccessful attempts |
| `software/` | Isolated h5py installation used to read legacy HDF5 files |

The previous dataset is `/data/pyoung3/supernova_data/50Am/analysis/50Am_plotting.out`; the new dataset is `/scratch/rpraka15/50Am_plotting_all_isotopes.out`. The filename `50AmPlotting.out` was not located. See `provenance/all_group_dataset_paths.csv` for alternate legacy paths. `final_chunks.out` is 1,438,981,268 bytes, or 1.44 decimal GB.

## Reproduction

Use `/home/rpraka15/.conda/envs/ccsn/bin/python`. Versions are recorded in `provenance/environment.json`; h5py 3.16.0 is isolated in `software/`. The existing SDF/time-series reader is `/scratch/rpraka15/supernova-grain-analysis/code/readers.py`; its source snapshot is retained under provenance. All expensive processing is submitted to Sol Slurm, account `grp_pyoung3`, partition `htc`, qos `public`.

Submit stages in dependency order, waiting for success and inspecting their logs:

1. `sbatch scripts/population.sbatch` parses, hashes and compares every plotting-file row.
2. `sbatch scripts/chunks.sbatch`, `sbatch scripts/raw_network.sbatch`, `sbatch scripts/provenance.sbatch`, and `sbatch scripts/group_inventory.sbatch` audit chunks, legacy HDF5 and filesystem provenance. The full group inventory enumerates about 21 million entries and is intentionally not a routine rerun.
3. `sbatch scripts/archive_probe.sbatch`, `sbatch scripts/eoc.sbatch`, `sbatch scripts/trajectories.sbatch`, and `sbatch scripts/initial.sbatch` inspect raw archive prefixes, every compact final-output CSV, every candidate trajectory and every archived initial composition. The trajectory stage requires the raw-network stage.
4. `sbatch scripts/figures.sbatch` and `sbatch scripts/synthesis.sbatch` use completed population/chunk/EOC results. Synthesis also requires the group inventory and archive probes.
5. `sbatch scripts/final_validation.sbatch` and `sbatch scripts/physical_source.sbatch` complete source-history, initial-composition and additional provenance checks.
6. `sbatch scripts/seeds.sbatch` compares every archived seed against the native initial composition after the physical-source stage succeeds.
7. `sbatch scripts/reports.sbatch` builds reports and the tracer-stage coverage ledger after all preceding stages succeed.

The audit caches are not general-purpose automatic invalidation caches. Original dataset receipts and end-of-audit stability checks identify the exact inputs used. For a new input revision, use a fresh audit directory and update input paths; do not silently reuse old caches. `chunk_audit.py` retries unsuccessful chunks but reuses successful receipts. No script runs a nuclear reaction network. Source scripts that move or delete files were inspected, never executed.

Run the small scientific failure-mode tests with:

```bash
cd /scratch/rpraka15/50AM_full_nucleosynthesis_comparison
/home/rpraka15/.conda/envs/ccsn/bin/python -m unittest discover -s scripts -p 'test_*.py' -v
```

## Numerical methods

Join unique, valid particle IDs; preserve original row order in parsed tables and sorted matched IDs in `alignment.npz`. Duplicate IDs are excluded from ambiguous joins. Missing columns, explicit `NA`, malformed tokens and numerical zero remain distinct.

Abundance equality: exact float equality first; then `abs(new-old) <= 1e-5*abs(old)` with zero absolute tolerance. Copied physical fields use relative tolerance `1e-9`. Raw final EOC versus six-significant-digit chunk values uses `5.1e-6`; B-weighted chunk values use `2e-5` to accommodate both rounded factors. Isotope number-ratio comparisons use `2e-5`. All numerical tolerances are computational comparison tolerances, not estimates of astrophysical uncertainty.

Relative differences and new/old ratios require nonzero old values; logarithms require both values positive. No abundance floor, pseudocount, clipping or forced normalization is introduced. Log differences are computed as `log10(new)-log10(old)` to avoid overflow of an intermediate ratio. Undefined values are NaN; overflow is counted separately.

Yields use particle mass in grams divided by `1.9889e33` g per solar mass. They are sums over explicitly stated available/matched cohorts; no unbound-ejecta selection, fallback correction or decay correction is silently applied. Available-isotope sums exclude elemental totals and count each physical nuclide once. See [the database guide](tables/README.md).
