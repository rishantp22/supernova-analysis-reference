# PRISM formats and prior CSV audit

Audited 2026-09-25. Findings describe inspected samples and current source, not a census or proof of the executable used by every historical job. No collaborator files were changed.

## One particle directory

A numeric directory such as `/data/pyoung3/groupdata/svierlin/jet3b_data/jet3b_outputs/10076` contains a postprocessed reaction-network calculation. Use `(model, particle_id)` as the key, never ID alone across models. Configuration paths identify inputs but may have moved since execution.

| File | Structure | Meaning |
|---|---|---|
| `control.json` | Inputs, network, output settings | Trajectory/composition paths, numerical start/stop, optional thresholds |
| `x.txt` | `Z A X` | Final isotope mass fraction, dimensionless |
| `y.txt` | `Z A Y` | Final abundance, `Y=X/A`; not particle mass or a normalized elemental number fraction |
| `xz.txt` | `Z X_element` | Final mass fraction summed over every network species with this Z |
| `yz.txt` | `Z Y_element` | Final Y summed over species with this Z |
| `xa.txt`, `ya.txt` | `A summed_value` | Sums over mass number A, mixing elements/isobars |
| `xtime.txt` | Repeated timestep/time/temperature/density header, followed by sparse `Z A X` rows | Composition evolution at accepted network output steps |
| `ytime.txt` | Corresponding time blocks and `Z A Y` | Abundance evolution |

Z is proton number, A total nucleon number, N=A-Z. `/scratch/svierlin/prism-v02/src/progression.f90:173` reads initial Z,A,X, computes N and assigns Y=X/A; line 199 normalizes the initial mass fraction. `src/output.f90:264` computes X=Y*A; line 332 sums all species for elemental X. The README at line 76 describes xZ as “Y(Z) mass fractions”; the implementation shows this is a typo.

The xtime timestep is an integration counter, not an SPH snapshot number or seconds. Header units are seconds, GK (multiply by 1e9 for kelvin), and g/cm^3. Parse header labels carefully because a label and number may touch. Never join these records to SPH files by row position.

## Sparse output is not zero

Defaults for all final/time X/Y streams are `lower_bound=1e-15` (`src/control.f90:200` through 207). JSON can override each stream (`src/datasets.f90:354`); the three inspected controls do not. Rows are emitted only for values **strictly greater** than the threshold (`src/output.f90:232`, 265, 335, 486).

An omitted isotope is not a measured zero. X and Y apply their cutoffs independently, so their species sets can differ. Missing rows can also mean outside the network, incomplete output, missing file, or an unverified association. Retain these statuses separately. Infer “at or below output threshold” only after verifying membership and successful relevant output; nonnegativity is an additional physical/validation assumption.

For final elemental totals prefer xz.txt: it sums **all network species before thresholding the total**. Summing only requested isotopes gives a subtotal. Summing sparse x.txt omits small species; an intermediate xtime sum also needs this omitted-tail qualification.

`/scratch/svierlin/prism-v02/input/extent/default.dat` uses **Z, N_min, N_max**, inclusive neutron-number bounds, not mass-number bounds (`src/networkSize.f90:148`). Example `6 2 19` permits carbon A=8 through 25. Initial species outside the extent may be mapped to a boundary (`src/progression.f90:181`), so initial input species and actual network membership are different concepts.

## Completion and time provenance

Final outputs are written only for `EVO_FINISHED` (`src/output.f90:230`, 262, 329), which occurs at `t >= t_final` (`src/progression.f90:391`). Time files may exist even if a run later fails (`src/output.f90:479`). Files may also be stale from an earlier execution. Preserve configured numerical stop, last saved time, final-file availability and logs/status evidence. Check final X against the last time block within printed precision, where that block actually reaches the endpoint.

The default `printout_rate=1` is at `src/datasets.f90:135`; other rates can leave the last saved block before the final state. x.txt itself has no timestamp. Do not treat all “last” rows as one common epoch or as fully decayed stable abundances.

All sample suffixes below are under `/data/pyoung3/groupdata/svierlin/`:

| Sample suffix | Config start s | Config stop s | Last xtime s | Blocks | X rows | Sum X |
|---|---:|---:|---:|---:|---:|---:|
| `jet3b_data/jet3b_outputs/10076` | 0 | 8000 | 8000 | 852 | 140 | 1.0000000104771447 |
| `cco2_data/cco2_outputs/118326` | 0.014655 | 860.69 | 860.69 | 65 | 151 | 0.9999999999999987 |
| `vconv_data/vconv_outputs/251189` | 0 | 1663.1 | 1663.1 | 71 | 148 | 0.9999999999999991 |

Their last headers are respectively xtime.txt lines 203397, 10439, 11268. These are sample checks, not completeness estimates. jet3b has 122 Y rows, with 18 X-only species; on common keys max absolute X-A*Y difference is 4.44e-16. Its sum(xz) is 1.000000010477149. Preserve that small normalization deviation in diagnostics rather than silently renormalizing.

All three configs reference `50Am_abund` initial compositions and model-specific trajectories. Shared initialization does not make the runs interchangeable. The jet3b control numerical stop at line 280 is 8000, while its comment says “10.0sec”; the other samples carry the same stale comment. Numerical settings and actual times take precedence.

## Previous scratch CSVs

Observed files are under `/scratch/rpraka15/final_chunks/`:

| File | Observed header | Interpretation |
|---|---|---|
| `chunk_1.csv` | `ID,category,actual_timestep,actual_time,A,B,C_sum` | Selected isotopes at three selected epochs; column A is actually Z, column B nuclear A |
| `final_chunks_B/chunk_1_B.csv` | `ID,category,A,B,mass_times_Csum` | Mass weighting; actual time/step dropped |
| `chunks_together_RP.csv` | `category,A,B,total_mass_times_Csum` | Aggregate across particles; IDs dropped |

`chunk_1.csv:2` is ID 100000, closest_to_100, step 77, actual time **122.085 s**, Z=1, A=1, zero. “Closest to 100” is not a measurement at 100 seconds.

Current collaborator scripts explain the schemas but do not prove byte-for-byte provenance of copied products. Paths below are under `/scratch/svierlin/prism-v02/`:

- `sum_eoc_manifest_RP.sh:14` uses eoc_manifest and historical all_outputs with no model column. Lines 78–90 select nearest saved times to 100/1000 seconds; lines 114–116 select the last. Lines 126–129 **replace missing selected-isotope rows with zero**, losing missingness. Line 100 adds rows by numeric time and isotope, so repeated blocks at the same formatted time would be summed. Line 131 uses `%g`, normally six significant digits.
- `sum_eoc_manifest_part2.sh:49` looks up mass by ID and numeric category, maps last to category 2000 at line 73, and substitutes zero for missing mass at line 77 (with a warning). Category 2000 does not prove abundance time is 2000 seconds. Actual time/step disappear at line 63. Its broad CSV glob at line 23 can also include aggregate CSVs if present.
- `sum_eoc_manifest_part3_RP.sh:239` sums mass-weighted values across IDs, then applies `1e-6` at line 251. Verify original mass units before calling totals solar masses or grams. The aggregate cannot provide per-particle data.
- `EOC_compile.sh:50` writes the clearer final schema Z,A,x,y and line 108 retains absent values as NaN. It still needs model/run/time provenance; lines 19–20 refer to historical all_outputs and 50Am_abund.

Use these products as leads and comparisons, not the sole authoritative particle dataset. Do not overwrite them or reinterpret inserted zeros as physical measurements.

## Prior Noh reports are not the actual 50Am explosion audit

`/scratch/rpraka15/output_data/DATA_STRUCTURE_ANALYSIS.md:5` explicitly describes a 1000-particle Noh test; line 12 names noh1k.0. `50AM_TEST_ANALYSIS.md:177` says 50AM_test was extracted from snsph.tgz, a code/test bundle. `ISOTOPIC_DATA_ANALYSIS.md:9` states that test lacks composition tracking. These conclusions must not be generalized to actual 50Am explosion outputs, their particle counts, layouts or composition availability.

## Remaining checks

Historical build identity, whole-run completion, model provenance of old archives, and mass units of old aggregates are unresolved here. This check read small source files and three particle examples; it did not extract full archives, scan million-directory trees, or submit jobs.

## Combined `.all.vrout` files contain repeated particle epochs

A bounded prefix scan of `/data/pyoung3/groupdata/cco2.all.vrout` found ID 23 at lines 1 and 944808; the two full row strings exactly match the first records of `/scratch/svierlin/prism-v02/cco2_vrout_files/r3g_1M_cco2_sph.1000.vrout` and `.1100.vrout`, respectively. In `/data/pyoung3/groupdata/vconv.all.vrout`, ID 0 repeats at lines 1 and 945629, exactly matching the first records of `/scratch/svierlin/prism-v02/vconv_vrout_files/run3g_50_1milvc2_sph.0000.vrout` and `.0020.vrout`. Each scan stopped at the first repeated ID. This supports concatenated snapshots at the beginning of each combined file; it does not establish the complete file order, absence of omitted snapshots, or all later block boundaries. No construction command was located in the bounded top-level script search.

`/data/pyoung3/gsvance/data_mining/cco2-SDF-reader.c:128` writes the headerless 14-column order **ID, x, y, z, temp, rho, pressure, vx, vy, vz, mass, Y_el, radius, speed**. Its variable named `vr` is the magnitude `sqrt(vx²+vy²+vz²)` at line 125, **not signed radial velocity**. These rows carry no time or snapshot identifier. Do not deduplicate on ID or join the combined file as a single-epoch particle table. Recover epochs from verified per-snapshot filenames and corresponding SDF `tpos` and unit factors, retaining snapshot identity in the row key; do not assume a constant particle count when splitting the aggregate.

The current trajectory builders use individual snapshot vrout files and external hardcoded time arrays: `generate_trajectories_cco2.sh:270` loops over 1000–19900 by 100, while `generate_trajectories_vconv.sh:250` loops over 0000–3380 by 20. Their assignments at lines 279/259 and output formatting at 280/260 are independent of the headerless combined files. These tables are provenance evidence, not authoritative substitutes for checking the original SDF time units.

The reader declares `Y_el` with a “for alignment” comment (`cco2-SDF-reader.c:66`, `SDF-reader.c:59`). This field's name and position alone do not establish that every source simulation initialized/evolved it as physical electron fraction. Keep the stored value explicitly unvalidated until its writer and consistency with the isotope-derived sum `Σ(Z/A)X` at the same epoch are established.
