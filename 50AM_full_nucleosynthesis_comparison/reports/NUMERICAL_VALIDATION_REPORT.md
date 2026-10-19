# Numerical and physical validation

## Population and physical quantities

The plotting comparison includes all 937,227 shared unique valid IDs. Every one of the 15 shared physical columns is exactly equal after numeric parsing, with no merely approximate or significant physical mismatch. This equality is constructed by the new exporter. An independent source check also confirms all 945,628 old-file IDs' final x/y/z, velocity, acceleration, mass, h and density against native SDF fields after conversion, at rtol=1e-5 and atol=0. Initial reported Ye also matches the native stored `Y_el` field at this precision.

Units: positions and h in cm; velocities cm/s; accelerations cm/s²; mass g; density g/cm³; peak temperature K; Ye dimensionless; isotope/element X dimensionless mass fractions. Native length, time and mass factors are 6.955e10 cm, 100 s and 1.9889e27 g. Legacy 'peak density' is density at the sampled peak-temperature point, not an independently maximized density. The old early-snapshot peak window differs from a maximum over all 61 snapshots.

**Correction to an early interpretation:** although the legacy exporter has a missing-record-to-zero path, the observed copied zero Ye values match zeros already in the original hydro `Y_el` field. Thus their presence is not attributed to an exporter lookup failure. There are 495,014 such zeros among valid new IDs. Native composition-derived Ye exposes their inconsistency; neither these zeros nor partial plotting Ye should be used to classify neutron-rich freeze-out.

All 3,052,610 legacy history times match physical snapshot times at rtol=1e-6, atol=1e-9 s. Temperature and converted density both match at 3,052,368 points with rtol=1e-10, atol=0. The remaining 242 points belong to 22 particles; their stored values include extreme signed magnitudes and one nonfinite point. Full particle flags and examples are retained, without repairing the HDF5 data.

## Mass closure and electron fraction

Available-isotope plotting sums exclude elemental totals and count each nuclide once. The inventories are incomplete, so deficits are not automatically errors. Over-unity sums are inconsistent with valid mass fractions even for incomplete inventories.

| dataset | particles | isotope_columns | rows_missing | over_unity_1e5 | sum_min | sum_p50 | sum_max | max_sum_id |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| old | 945628 | 65 | 0 | 279 | 0 | 0 | 1.374807 | 157594 |
| new | 937227 | 223 | 937227 | 1 | 2.471662e-05 | 0.9999542 | 1.977369 | 7753 |


All 524 legacy raw species give sum-X range 3.352899894e-07–1.38239028; 2251 particles exceed 1.00001 and 2757 fall below 0.99999. Maximum is ID 157594. These are checks of saved results, not a successful-solver certificate.

The full independent EOC archive gives finite-X sums from 0.99999489274 to 1.00000790957; all 938,732 are within 1e-5 of unity. Its composition-derived Ye ranges 0.4811954801–0.4995180553. Censored/missing X remains recorded. Closure cannot validate trajectory units, seeds, reaction rates or historical run identity.

Native initial composition uses slots 1–20; the source explicitly excludes duplicate free-nucleon slots 21–22. Native composition closes within approximately 8e-6, while native stored Ye disagrees severely in much of the outer population. Archived PRISM seeds have sums 0.999232410975–1.000814473939 and are normalized explicitly by the network initialization. The audit preserves original seed values and separately labels normalized diagnostics.

All new elemental totals (59 columns) are unavailable. Legacy elemental totals are compared with available isotope subtotals in `elemental_total_checks.csv`; a subtotal deficit can reflect unexported isotopes or thresholded queries and is not sufficient evidence of incorrect elemental mass. No new total is synthesized and substituted into the original file.

## Ratios, decay and censoring

Number ratios use `(X_i/A_i)/(X_j/A_j)`. The 172 audited pairs include every available isotope against a declared within-element reference plus all requested Mg–Al, Ca, Ti, Cr and Fe pairs. Missing columns and zero denominators remain explicit; no pseudocount or abundance floor is used. The all-particle ratio files and quantiles are in `tables/ratios/` and `number_ratio_statistics.csv`.

The old/new tables are not demonstrated to be at one common freeze-out or decay epoch. Chain subtotal checks test whether *isolated* decay alone could explain a difference; they do not apply decay corrections and do not prove decay when a sum agrees. The 44Ti→44Sc→44Ca chain is unavailable as a complete plotting chain because 44Sc is omitted. The 26Al output has no state separation. No inferred initial Solar System 26Al/27Al is equated with a model mass fraction or model number ratio.

| chain | old_missing | new_missing | valid_closed_chain_pairs | inconsistent_with_isolated_decay_1e5 |
| --- | --- | --- | --- | --- |
| 26Al -> 26Mg | unavailable | unavailable | 937227 | 935935 |
| 41Ca -> 41K | unavailable | unavailable | 937227 | 937226 |
| 44Ti -> 44Sc -> 44Ca | 44Sc | 44Sc | 0 | 0 |
| 53Mn -> 53Cr | 53Mn | unavailable | 0 | 0 |
| 56Ni -> 56Co -> 56Fe | 56Co | unavailable | 0 | 0 |
| 57Ni -> 57Co -> 57Fe | 57Ni;57Co | unavailable | 0 | 0 |


Thresholded source output is different from exact zero. The EOC comparison records 135,360,101 plotting-zero cells without a finite final-X counterpart and 5,682 plotting-positive cells without one. These cases cannot be forced into ordinary relative/log comparisons. Source X=A Y consistency holds for every finite pair tested.

## Tolerances and numerical limits

The primary abundance comparison uses rtol=1e-5, atol=0; copied physical fields rtol=1e-9, atol=0; isotope number ratios rtol=2e-5, atol=0. Exact matches are separately counted. Raw query rounding checks use 5.1e-7; EOC/chunk output rounding 5.1e-6; B weighting 2e-5. All are computational tolerances, not physical error bars. The candidate density export uses rtol=5.1e-6 plus atol=5.1e-6 native units to accommodate five printed decimal places; without that documented rounding allowance, small native densities appear discrepant.

All statistics use the full relevant population. Distribution percentiles omit undefined/nonfinite values and accompany availability/overflow counts. Log differences are evaluated only for positive pairs without clipping tiny abundances. No input normalization is imposed except the explicitly labeled diagnostic reproducing the network's own initial normalization.
