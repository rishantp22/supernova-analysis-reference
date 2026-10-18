AM chunk export: final_chunks.out

Rows: 937,229; columns: 298; delimiter: TAB; missing value: NA.
One header row, exactly the order/spelling in List of isotopes_elements_variables.txt.
One row per ID represented in the supplied chunk files; no invented IDs.
Rows are ordered by numeric chunk number, then numeric ID within each chunk.
Id is the source simulation/SPH particle ID (the user's grain/clump identifier).

Composition category: last. Isotope columns copy C_sum as dimensionless
mass-fraction values from final_chunks/chunk_N.csv, without normalization or decay.
The two chunk sets were joined/checked using ID, category, proton number (source A),
and mass number (source B). B contains particle-mass-weighted isotope values;
these are checked separately and are NOT substituted for isotope mass fractions.
Mass comparison tolerance is 2e-5 relative because source products are rounded.
B mass discrepancies: 0 of 208064838 checked cells.
B values are in native code mass, NOT already solar masses: multiply by 1e-6
using the convention documented in local SNSPH code and the original aggregation
script. This factor is source-code-supported, not declared in the AM SDF header.
mass_fraction_validation.json gives totals for each category separately, counts
of per-ID isotope sums close to one, and tests of both possible B unit hypotheses.
normalization_by_id.csv contains the sums for every ID and all three categories.
These sums cover only the selected isotopes and supplied IDs; they do not by
themselves establish total ejected mass. Never add the three epochs together.

Physical fields come from run3g_50Am6_sph.2000.out, joined by ID. All matching
text values were compared with the raw AM .2000 SDF to printed precision.
Position, velocity, acceleration, mass, h, Density retain original CODE UNITS.
Y_{e} retains the stored Y_e value; this is not a recomputed/validated electron
fraction. No CGS conversion is applied. Hydro iteration is 2000 and tpos is
1575.58 code time. Composition actual_time is retained separately in grain_epochs.csv;
matching IDs does not imply matching physical epochs. Historical trajectory time
and density problems described in the existing audit are not repaired by reshaping.

Peak temp is the maximum stored temp (K) and Peak density the maximum stored rho
(code density) across 61 saved raw AM snapshots. These are separate
maxima, not necessarily simultaneous, and are sampled hydro maxima, not guaranteed
continuous peaks or maxima of the PRISM input trajectories. Per-ID snapshot counts
are in grain_epochs.csv; complete source coverage is in peak_sources.json.

Element policy: na. With 'na', complete elemental totals are unknown
because chunks contain only selected isotopes. With 'subtotal', element columns
sum only the available selected isotopes and must be interpreted as subtotals.
An element with no selected isotope remains NA. X_{n} is absent from the chunks
and remains NA; hydro-network neutron fractions are not mixed into PRISM composition.
Columns entirely unavailable: X_{n}, H, He, Li, Be, B, C, N, O, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ge, Se, Kr, Sr, Y, Zr, Nb, Mo, Ru, Pd, Ag, Cd, In, Sn, Te, I, Xe, Ba, La, Nd, Sm, Eu, Gd, Dy, Er, Hf, Ta, W, Os, Pt.

Existing zero C_sum values are preserved exactly as numerical zeros. The historical
chunk generator inserted zeros for absent isotope rows, so these zeros do not
establish a measured physical zero. No new missing values are filled with zero.
The chunk copy script points to historical collaborator products without an explicit
model column. The user's AM attribution and physical ID/mass crosschecks are recorded;
this export is not a new validation of every original PRISM run's provenance.

validation.json contains counts, missingness and final SHA-256. chunk_validation.json
contains all per-chunk checks and SHA-256 fingerprints of both source chunk sets.
grain_epochs.csv preserves per-ID composition times, hydro coverage, isotope sums,
and counts of source zeros. Output precision for parsed abundance/peak values is
10 significant digits; original physical text values are retained as written.
All source files were opened read-only and left unchanged. No source dataset was deleted.
