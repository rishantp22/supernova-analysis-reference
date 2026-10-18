50Am_plotting_all_isotopes.out

937,229 unique particle IDs; 298 columns; one header line; comma-space separated.
Rows are sorted by numeric ID. The original requested list determines column order.
Names follow the example: id, density, peak temp, peak density, X_{H}, X_{1H}, etc.
column_map.csv maps every original requested name to the output name and source.

The first 16 physical fields are copied by ID from the actual reference file:
/data/pyoung3/supernova_data/50Am/analysis/50Am_plotting.out
Positions and h: cm; velocities: cm/s; acceleration: cm/s^2; mass: g;
density and peak density: g/cm^3; peak temp: K. Y_{e} is the PROGENITOR
electron fraction used by the reference generator, not final snapshot Y_e.
Peaks follow the reference's early-snapshot definition and printed precision.
The length/mass/density conversions agree with the validated scratch physical
file for all matched IDs. Local SNSPH conversion conventions are preserved.

All 222 requested isotope columns supplied by the chunks retain their last-category
C_sum values. No isotope value was replaced by a legacy burn-query abundance.
The reference plotting file is an older burn-query product; its composition is
different from the supplied PRISM chunks. Mixing its elemental totals with the
new isotope values would not form one consistent composition.
Consequently all 59 elemental-total columns and X_{n} remain NA, as their
complete values cannot be recovered from the selected-isotope chunks. NA is
missing information, never zero. Existing source zeros are preserved.

The reference has 945,628 IDs, of which 8,399 have no chunk abundances.
The output contains every one of the 937,229 chunk IDs. Example IDs 0, 2 and 3
are among those missing from the chunks, so they cannot be extended with the
requested chunk isotopes. The prior hydro_ids_without_chunks.csv lists those IDs.

The isotope sum anomaly at ID 7753 remains unchanged. No renormalization or decay
correction is applied. Abundance epochs are retained in the prior grain_epochs.csv.
All source files, the original plotting table and final_chunks.out remain unchanged.
