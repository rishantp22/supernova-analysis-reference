50Am_plotting_complete.out

All 945,628 AM IDs from the actual 50Am_plotting.out reference are included,
sorted numerically, with 298 comma-space-separated columns. Header spelling follows
the example (id, density, X_{4He}, etc.); column order follows the requested list.

The first 16 fields are copied by ID from the reference plotting table, in CGS units.
Y_{e} is its progenitor electron fraction. Peaks retain its early-snapshot definition.
937,229 IDs have all 222 requested isotope values supplied by the last-category
chunks. Those values are unchanged from the validated chunk export.
8,399 IDs, including 0, 2 and 3, have no supplied chunk compositions: their
physical fields are present and all 282 composition fields are explicitly NA.
complete_model_coverage.csv identifies these rows.

For all rows, 59 full elemental totals and X_{n} remain NA because the selected
isotope chunks cannot determine them. NA is unavailable information, not zero.
The older reference's composition belongs to a different abundance calculation;
it was not mixed into the new chunk isotope values. No missing data were invented.
All 298 requested variables are represented, but a fully populated numerical table
cannot be obtained from these chunks alone.

Original sources and earlier outputs remain unchanged. See column_map.csv for the
field mapping and README.txt for physical definitions and source distinctions.
complete_model_validation.json records counts, checks and SHA-256.
