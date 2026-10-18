# Results of the 25 September 2026 structural audit

The audit maps all **298 requested fields** across four observed model labels,
giving 1,192 rows in the source/field map. It inspected **93 PRISM particle
examples**: six 50Am archive examples, 27 jet3b, 28 cco2 and 32 vconv. These are
deterministic structural examples, not a statistically representative sample.
The requested-value sample has 27,714 rows including explicit missing statuses.

All 483 selected-directory SDF headers parsed and matched their declared record
count, schema and file length. The selected final snapshots contain:

| Observed model | Iteration | Particle records | Duplicate IDs |
|---|---:|---:|---:|
| 50Am | 2000 | 945,628 | 0 |
| jet3b | 8000 | 945,628 | 0 |
| cco2 | 19900 | 931,564 | 0 |
| vconv | 4200 | 945,342 | 0 |

The differing counts, populations and epochs must remain visible when joining
data. The `vconv` association with the requested vconvL is not fully established.
See `tables/snapshot_catalog.csv` and `tables/hydro_diagnostics.csv` in scratch.

## Read gaps are findings

The broader sample found six jet3b entries without `xz.txt`: IDs 10354, 103110,
103295, 103480, 103664 and 103849. Their elemental values are marked missing
source, not zero. Isotope data cannot recover an exactly equivalent elemental
total from sparse x.txt without acknowledging the omitted contribution.

cco2 particle 14 lacks the inspected x, xz, y and xtime files; its configured
trajectory contains no data rows. Its normalization and electron fraction are
blank, not zero, and no endpoint validation is claimed. The vconv particle 14
xtime file exceeded the 32 MiB per-file sample cap; its endpoint is unverified,
which is not evidence that the job failed. Six configured 50Am trajectory paths
are missing; alternative archived candidates were inspected separately without
assuming they are the actual historical inputs.

The run therefore reports `completed_with_read_errors`: 11 sample-file read/cap
records and seven input-trajectory error records, with no snapshot-header
errors. Exact paths and reasons are in `debug/sample_read_errors.json` and
`debug/prism_errors.json`. These records distinguish unavailable data from
parser failures and intentional bounds.

In the inspected outputs with a readable saved history, all checked endpoints
match their configured stop within printed precision. This is an endpoint
check, not validation of the physical trajectory or job/build provenance.
The largest inspected final PRISM X sum is about 1.000011325 (jet3b ID103);
small deviations are preserved instead of silently renormalized.

## Most consequential findings

- **Input units:** observed trajectory values corroborate the time/density
  discrepancies in current scripts. Existing cco2/vconv times lack the header's
  ×100 conversion; jet3b uses iteration as time. Their density factor is about
  4.163 times smaller than the reference conversion. 50Am archived candidates
  have additional temperature/history issues. See `reports/unit_crosscheck.md`
  and its numerical CSV in scratch; 50Am/jet3b reference factors are inferred.
- **Archive completeness:** 50Am batches 1, 8 and 9 have matching no-space
  creation failures in logs. Full integrity is unverified. Prefix reads cannot
  certify the remainder of a compressed archive.
- **Stored field quality:** 20,000-record late jet3b/vconv samples have zero
  stored Y_el. Sampled 50Am/jet3b hydro composition-slot sums reach roughly 1.03.
  Interpretation requires run/network provenance, not an automatic repair.
- **Prior products:** bundled Noh-test reports and old mass-aggregated CSVs do
  not supply the requested authoritative per-particle explosion table.

## Delivered artifacts and next scientific decision

The scratch report `reports/index.html` includes seven diagnostic figures,
searchable field mappings, sampled values and links to inventories. Code,
documentation and reader tests are in the local Git repository; code snapshots
and hashes are also retained in scratch. Source data were read-only throughout.

A full all-particle science export remains a separate stage: establish the
correct historical inputs or rerun affected postprocessing, resolve the vconvL
association, choose abundance epochs/decay treatment, and establish completed
particle coverage and archive integrity. These are concrete open issues found
by this audit, not assumptions hidden in a table of apparently complete values.
