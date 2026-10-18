Read-only comparison of final_chunks.out and 50Am_plotting_all_isotopes.out

Original paths:
/scratch/rpraka15/final_chunks.out
/scratch/rpraka15/50Am_plotting_all_isotopes.out

compare_outputs.py waits for at least the prior complete first-file size and 30 seconds
of unchanged size/mtime/ctime/inode before starting. It checks cmp first, then
hashes and validates every input row while writing temporary records. It uses
1024-row NumPy blocks and GNU sort with a 128 MiB buffer, one thread, and scratch
temporary files. The original files are opened read-only. Metadata is checked
again after each scan and at the end; changes invalidate the comparison.

Rows match by integer particle ID. Headers map by name: Id/id, Density/density,
Peak temp/peak temp, Peak density/peak density, and isotope/element names such as
56Fe/X_{56Fe}. Y_{e} retains its name; X_{n} maps consistently to n. Duplicate IDs
are counted and shared duplicate groups are excluded from ambiguous comparisons.
Missing IDs, malformed records, column counts, ordering, missing cells, zero and
negative counts are recorded. NA is missing information, never a numeric zero.

Every numeric pair is checked for exact float64 equality before conversion.
All composition token strings are also checked for exact equality after delimiter
normalization when the mapped column order is the same. Unit conversions are
taken from the local export source make_plotting_format.py, not fitted to data.
After conversion, abs(final * factor - plotting) <= 1e-300 + rtol*abs(plotting)
is checked at rtol=1e-9 and rtol=1e-5. The latter separates the source's six-digit
printing from larger semantic differences. Exact isotope differences remain
reported independently, including subnormal values; tolerances never replace NA.

The saved comparison.json, column_comparison.csv, mismatch_examples.json and
selected_particle_examples.json contain the measured result. Column counts are
delimiter-aware: whitespace alone incorrectly splits the two peak-field names.

Successful Slurm comparison job: 66028953 (one CPU, 2 GiB requested, 1-hour limit).
Job 66028845 was cancelled while waiting: the stable input is two bytes larger
than the prior manifest. The final job audits that discrepancy rather than
requiring byte identity with the prior manifest before performing the audit.
Job 66028903 stopped when it encountered two invalid numeric tokens containing y.
The final comparator reports such invalid cells separately from NA and continues
comparing the remaining valid cells. Malformed IDs are reported and excluded from
integer-ID matching without attempting repairs.
Submission script: compare_outputs.sbatch. Test fixtures: test_compare_outputs.py.
Prior export manifests supply expected sizes/hashes; this job independently
recomputes those hashes and all pairwise data checks.

The initial 19,660,801-byte observation was during an ongoing file transfer;
it is not the complete size. A subsequent size decrease indicated a restarted
transfer. Only a stable file at least as long as the prior complete export is
accepted. A second checksum excluding blank lines diagnoses added blank lines
without changing the actual file hash or suppressing malformed-row reporting.
