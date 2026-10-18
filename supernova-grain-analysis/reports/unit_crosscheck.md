# Existing trajectory unit crosscheck

This check reads existing trajectory files referenced by audited particle configurations. Particle IDs are selected explicitly from the raw SDF records. Epoch matches require temperature agreement within text-rounding tolerance and agreement with the inspected generation script time convention. Rows are never joined by position. Some early vconv temperatures repeat, so temperature alone would be insufficient.

| Model / variant | Particle | Iteration | Written time (s label) | Reference seconds | Density ratio |
|---|---:|---:|---:|---:|---:|
| jet3b / existing_file | 10076 | 0 | 0 | 0 | 0.2401969 |
| jet3b / existing_file | 10076 | 100 | 100 | 0.217 | 0.2401971 |
| jet3b / existing_file | 10076 | 8000 | 8000 | 9.64535e+07 | 0.2401968 |
| cco2 / existing_file | 118326 | 1000 | 0.0147 | 1.4655 | 0.2401964 |
| cco2 / existing_file | 118326 | 1100 | 0.0167 | 1.6675 | 0.2401966 |
| cco2 / existing_file | 118326 | 19900 | 861 | 86068.8 | 0.2401964 |
| vconv / existing_file | 251189 | 0 | 0 | 0 | 0.2401962 |
| vconv / existing_file | 251189 | 20 | 3.5e-06 | 0.00035 | 0.2401963 |
| vconv / existing_file | 251189 | 3380 | 1660 | 166313 | 0.240197 |
| 50Am / archived_original | 191136 | 0 | 0 | 0 | 169152.6 |
| 50Am / archived_original | 191136 | 20 | 20 | 0.02 | 169152.4 |
| 50Am / archived_original | 191136 | 2000 | 2000 | 157558 | 169238.6 |
| 50Am / archived_tempcorr | 191136 | 0 | 0 | 0 | 169152.6 |
| 50Am / archived_tempcorr | 191136 | 20 | 20 | 0.02 | 169152.4 |
| 50Am / archived_tempcorr | 191136 | 2000 | 2000 | 157558 | 169238.6 |

The sampled existing files corroborate the script-level discrepancies: jet3b writes iteration as time; cco2 and vconv write rounded code time, missing the factor 100; all three use approximately 1.42e-6 for density instead of 5.911824660348532e-6. The density ratio is approximately 0.2401965 (reference is 4.16326 times larger). These are observed output discrepancies, not only potential script hazards. They do not by themselves establish the downstream effect on nuclear abundances.

50Am particle 191136 has a missing configured trajectory path; two archived candidate variants were read directly from ZIP members without extraction. The original variant writes Kelvin values under a GK label (factor 1e9); the tempcorr variant fixes that temperature conversion. Both store iteration as seconds and raw code density under the CGS label. Archived variants have not been established as the actual input consumed by the sampled PRISM run. The jet3b reference factors are inferred from shared source code because the sampled headers omit conversion factors; the same limitation applies to 50Am. cco2/vconv factors are explicit in their headers. The vconv samples here do not resolve vconvL model identity.

Detailed numerical provenance: `tables/unit_crosscheck.csv`. Missing or ambiguous matches: `debug/unit_crosscheck_notes.json`.
