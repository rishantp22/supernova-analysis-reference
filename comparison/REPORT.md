# Comparison of the current scratch tables

The files are **not byte-for-byte identical**: `cmp` reports the first difference
at byte 1, line 1. They share identical isotope data for every valid shared ID,
but differ in physical units, some field definitions, ordering, and input errors.

| Property | `final_chunks.out` | `50Am_plotting_all_isotopes.out` |
|---|---:|---:|
| Full path | `/scratch/rpraka15/final_chunks.out` | `/scratch/rpraka15/50Am_plotting_all_isotopes.out` |
| Current bytes | 1,438,981,268 | 1,766,352,794 |
| Data rows, excluding header | 937,229 | 937,228 |
| Fields per row | 298 | 298 |
| Unique valid integer IDs | 937,229 | 937,227 |
| Duplicate valid IDs | 0 | 0 |
| Invalid IDs | 0 | 1 |
| Invalid non-ID numeric cells | 2 | 0 among valid-ID rows |
| Delimiter | Tab | Comma and space |

There are no blank lines, wrong-width records, or duplicate normalized column
names. Both headers have the same 298 fields in the same order after mapping
`Id` to `id`, physical-name capitalization, and isotope/element names such as
`56Fe` to `X_{56Fe}`. Names with spaces are parsed using their delimiters.
Final-chunk rows follow chunk ordering (937 numeric-ID inversions); plotting
rows with valid IDs are in ascending numeric order.

## Input errors and missing IDs

These values were observed and left unchanged:

| File | Line | ID | Field | Literal value |
|---|---:|---:|---|---|
| final_chunks.out | 16 | 10016 | vz | `-0.045y2421` |
| final_chunks.out | 30 | 10034 | ax | `7.02723e-0y5` |
| plotting file | 2 | invalid | id | `y193` |

IDs **23 and 193** exist in the final-chunk table but have no valid integer-ID
match in the plotting table. The literal `y193` was not silently changed to 193.
There are no valid IDs found only in the plotting file. All numeric comparisons
below use the **937,227 valid shared IDs**; invalid numeric cells are reported
separately and excluded from arithmetic, not counted as `NA` or zero.

Both files changed during the early inspection. The first file initially appeared
as 19.7 MB while a transfer was in progress; the final completed file is 1.44 GB.
Both current SHA-256 values differ from the saved export manifests. The two
extra bytes in final_chunks are not added blank lines: two numeric tokens contain
an extra `y`. The plotting file is now 1,572 bytes smaller than its saved manifest.
The final audit waited for stable metadata and verified it remained unchanged
throughout the full scan and comparison.

## Isotopes and missing values

All **222 isotope columns**, containing **208,064,394 compared numeric values**,
match exactly. All composition token strings also match exactly after delimiter
normalization for every shared valid ID, so this is not a tolerance-based result.
For example, ID 10000 has `16O = 1.97065e-15` and `56Ni = 0.123336` in both files.

The 59 elemental-total columns and neutron column (`X_{n}`) are `NA` throughout
the valid-ID records in both files. There are no missing-value disagreements on
matched IDs. These 60 unavailable columns are kept distinct from numeric zeros.

## Units and physical values

Headers do not declare units. The following factors and definitions are recorded
in the local exporter `make_plotting_format.py` and its source documentation.

| Field group | Multiply final-chunk values by | Plotting units |
|---|---:|---|
| x, y, z, h | 6.955e10 | cm |
| vx, vy, vz | 6.955e8 | cm/s |
| ax, ay, az | 6.955e6 | cm/s² |
| mass | 1.9889e27 | g |
| density, peak density | 1.9889e27 / (6.955e10)³ | g/cm³ |
| peak temp | 1 | K |
| Yₑ, isotope mass fractions | 1 | dimensionless |

All valid numeric pairs in the 12 ordinary physical columns agree after these
conversions at `rtol=1e-9`, `atol=1e-300`. For example, ID 10000 has raw
`x = 439.069`, which converts exactly to the plotting value `3.053724895e13 cm`.
The invalid `vz` and `ax` cells remain exceptions.

Three fields also differ in definition. Final-chunk peaks are independent maxima
over 61 saved hydro snapshots; plotting peaks follow the legacy reference's
early-snapshot definition. Final-chunk `Y_{e}` is the stored final hydro value,
whereas the plotting reference uses progenitor electron fraction.

| Field | Mismatches after conversion, rtol=1e-9 | Mismatches after conversion, rtol=1e-5 |
|---|---:|---:|
| Peak temperature | 933,037 | 536,333 |
| Peak density | 931,942 | 559,487 |
| Yₑ | 819,392 | 810,789 |

For both thresholds, the exact rule is
`abs(final * factor - plotting) <= 1e-300 + rtol * abs(plotting)`.
The second threshold accommodates six-significant-digit printing; values still
outside it are not explained by that printing precision. It does not prove that
every smaller difference is purely rounding. Exact isotope checks are separate.

Examples distinguish the effects:

- Rounding-scale peak temperature difference, ID 1167: `3,827,135,744 K` versus
  `3,827,140,000 K` (4,256 K; about 1.11e-6 relative).
- Larger peak-density difference, ID 1910: converted final-chunk value
  `3,361,604.350505526 g/cm³` versus `3,095,839.30806 g/cm³` (about 7.91% relative
  to the larger value), which remains after unit conversion.
- Electron fraction, ID 1167: `0.590934` versus `0.499970`.

## Audit record

Slurm job **66028953** completed with exit code 0 in **4 min 28 sec**, requesting
one CPU and 2 GiB. The Python process reported 71.3 MiB maximum RSS; Slurm's batch
step reported 153,384 KiB including the external sort. Sorting used a 128 MiB
buffer, and numeric calculations used 1,024-row blocks. Four small correctness
tests passed. Both inputs were opened read-only and their metadata was unchanged
during the successful audit.

Current SHA-256:

```text
final_chunks.out
bf6d52834b59d316c439f40d1e2bf36ab4f14d4fee7c535331deaca62e3bf290

50Am_plotting_all_isotopes.out
a5ad5ce7dc7a4beb8e75292fa0a4c64c56d71dde4bc16c8cedc223d554a55f24
```

Full measurements and examples are in [validation/comparison.json](validation/comparison.json),
[validation/column_comparison.csv](validation/column_comparison.csv), and the
other small validation files in that directory.
