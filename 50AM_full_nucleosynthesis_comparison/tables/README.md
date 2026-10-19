# Complete comparison database

All arrays are ordinary NumPy NPY/NPZ files, readable without a database server. Use `np.load(path, mmap_mode='r')` for NPY arrays and process one species at a time. NPZ files are compressed and decompress on access; each per-species file contains the full matched population, not a sample.

| Files | Row key and meaning |
| --- | --- |
| `old_values.npy`, `new_values.npy` | Original-file row order, float64 column-major arrays; columns and raw headers in corresponding `provenance/*_scan.json` |
| `old_ids.npy`, `new_ids.npy` | Integer particle IDs; invalid ID represented by −1 only as an internal sentinel; original malformed token retained in `*_invalid_cells.csv` |
| `alignment.npz` | Sorted shared unique valid `id`, with `old_row`, `new_row` indices |
| `species/NUCLIDE.npz` | Every aligned particle: `old`, `new`, `delta`, `relative_delta`, `ratio`, `log10_ratio`, `valid_pair`, `significant` |
| `species_inventory.csv` | Every plotting-file column, canonical nuclide, original name, availability and finite counts |
| `species_statistics.csv`, `physical_statistics.csv` | Counts, exact/approximate/significant matches, percentiles, extreme IDs and mass-weighted yields |
| `yields_all_available_particles.csv` | Separate populations, including unmatched valid rows; unavailable species yield is NaN |
| `yield_coverage_decomposition.csv` | Yield contributions from IDs with and without an available legacy raw result |
| `old_closure.npz`, `new_closure.npz` | Partial available-isotope sums, unaccounted mass, partial Ye, reported progenitor Ye and missing-column counts |
| `ratios/NUM_over_DEN.npz` | Number ratios on the same `alignment.id`; unavailable/zero-denominator cases remain NaN |
| `number_ratio_statistics.csv` | Complete 172-pair ratio audit; status records distinguish missing columns and zero denominators |
| `chunks/chunk_N.npz`, `all_chunk_particles.npz` | ID, three category times/steps, selected sums, partial Ye, zero counts and plotting correspondence |
| `raw_old/out50N/mass_fractions.npy` | Every stored legacy particle × all 524 nuclear species; row IDs and Z/A in `metadata.npz` |
| `raw_old/out50N/histories.npz` | Ragged historical time/temperature/density arrays: `id`, `offset`, `time_temp_rho`; units s, K, g/cm³ |
| `old_raw_vs_plotting.csv` | Direct raw-to-legacy-export comparison, including tiny positive raw values exported as zero |
| `old_raw_vs_new_selected_species.csv` | New selected species versus raw legacy results on 269,113 shared IDs |
| `eoc/block_NNNNN.npz` | Complete ragged EOC archive: ID, offsets, Z, A, X, Y; explicit NaN retained |
| `eoc_values.npy`, `eoc_ids.npy` | Dense memory-mapped EOC X view, columns in `eoc_species_inventory.csv` order |
| `eoc_row_present.npy` | Boolean matching dense EOC shape: false means absent isotope row; true with NaN X means explicit unavailable X |
| `old_raw_vs_new_eoc_all_shared_species.csv` | All 262 shared raw-network species compared on 270,706 shared raw/EOC IDs, finite availability per species stated |
| `eoc_plot_mismatch_cells.csv` | Every finite EOC/plotting mismatch with particle, isotope, both values and ratio |
| `trajectory_population.csv.gz` | All 943,900 candidate histories, source matches, physical and as-labeled cooling diagnostics |
| `initial_composition_population.csv.gz` | All 943,900 archived seed files, closure, composition Ye and explicitly labeled normalization diagnostic |
| `hydro_initial_composition.npz` | Native initial SDF species and fractions, with unused duplicate slots 21/22 excluded according to `unburned.c` |
| `prism_initial_X.npy`, `prism_initial_ids.npy` | Every archived PRISM seed mass fraction, without normalization; columns follow `initial_species_inventory.csv` |
| `prism_seed_vs_native_species.csv`, `prism_seed_vs_native_particles.csv.gz` | Full 943,900-particle native/normalized-seed comparisons; normalization is explicitly diagnostic |
| `tracer_stage_coverage.csv.gz` | Union of discovered tracer IDs with stage-availability flags; never infers missing results to be zero |

`all_new_species_environmental_bins.csv`, `thermodynamic_binned_statistics.csv`, and `all_species_physical_cooling_bins.csv` include all available species and the complete relevant populations. Environment plots based on old plotting fields are descriptive; copied zero peak/Ye metadata must not be interpreted as physical burning conditions. The archive-based physical cooling time is the first sampled fall below peak/e, not a fitted expansion timescale. No crossing produces NaN.

`eoc_species_inventory.csv` includes rows seeded into the EOC CSV inventory even when final X is unavailable. It is not the full configured PRISM network. A finite-value EOC yield is a sum of recorded contributions; omitted/censored states are not certified exact zeros. `old_raw_species_inventory.csv` describes the complete stored legacy network; `prism_default_extent_species.csv` describes the current configured extent. `complete_network_export_inventory.csv` combines network and export availability without filling missing species with zero.

Example: all-particle comparison of 26Al, with zero/missing cases preserved:

```python
from pathlib import Path
import numpy as np
root = Path('/scratch/rpraka15/50AM_full_nucleosynthesis_comparison')
ids = np.load(root/'tables/alignment.npz')['id']
with np.load(root/'tables/species/26Al.npz') as d:
    valid = d['valid_pair']
    changed = valid & d['significant']
    print(ids[changed], d['old'][changed], d['new'][changed])
```

The original-file scan and per-archive member hash tables provide provenance back to source bytes. No output array should be matched to another by row position unless their documented row key is identical.
