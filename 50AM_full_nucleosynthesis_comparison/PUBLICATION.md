# GitHub reference and full export

This folder records the completed 50AM audit from `/scratch/rpraka15/50AM_full_nucleosynthesis_comparison`.
The destination is the existing private repository
`rishantp22/supernova-analysis-reference`, in its own new
`50AM_full_nucleosynthesis_comparison/` directory. It does not change the presolar oxide repository.

The original reference includes 139 scientific code/report/table files.
All 226 generated figures were subsequently added at the user's request.
They include all audit Python/Slurm scripts, five reports, the database guide,
and top-level summary CSV/JSON records no larger than 2 MiB apiece.
All scripts, reports and tables are byte-identical to the completed audit;
the README adds a publication-scope notice. Both source and published hashes are recorded.
`PUBLICATION_MANIFEST.csv` records each original path, size and SHA-256.
`EXCLUDED_FILES.csv` lists every original file omitted from Git, with its size.
The publication metadata and this folder's `.gitignore` were added separately.

The complete computer export retains every file, including the full-population
databases, 113 PNG/PDF figure pairs, logs, source snapshots and local h5py copy.
It is prepared on Sol as
`/scratch/rpraka15/50AM_export_20261018/50AM_full_nucleosynthesis_comparison.tar.gz`.
Its `export_receipt.json` and SHA-256 sidecar record the final size and verification.
All 113 PNG and 113 PDF figures are included in GitHub (23,726,499 bytes).
No archive or full numerical database is uploaded to GitHub.
Original datasets elsewhere on Sol are not bundled in the audit folder.

## Reproduction and portability

Read `README.md` for stage order and `provenance/environment.json` for versions.
Scripts preserve the validated Sol paths; this is a research reference, not a
portable installed package. Full calculations require the original Sol inputs
and a writable scratch output directory. Slurm scripts retain that directory.
The shared SDF reader is also archived elsewhere in this repository at
`../supernova-grain-analysis/code/readers.py`.
The synthetic aggregation test reads the historical
`/scratch/svierlin/prism-v02/sum_eoc_manifest_RP.sh`; its snapshot is in the full export.
Source snapshots and nuclear controls remain available in that export.
Missing database paths in reports refer to the complete export, not to
files silently omitted from the numerical investigation.

Read `reports/KEY_FINDINGS.md` before scientific use: neither plotting dataset
has sufficient validation for quantitative publication. The audit does not
run a new nuclear network or correct the original abundances.

## Source changes observed during export preparation

The code/report/table commit `ac27975003cd8d67a6ca3feb628c43fb6c83e285`
preserves the audited results before subsequent scratch-file edits.
During archive preparation, `tables/eoc_plot_mismatch_particles.csv` changed:
literal `y` characters appeared at the ends of numeric cells for IDs 9460 and
10559. These are not scientific recalculations. The original scratch file was
left untouched; its edited contents are not substituted into the validated Git
table. The full export records the exact captured contents and any changes
detected during packaging. Consult its manifest and receipt when comparing
that export with this Git snapshot.
