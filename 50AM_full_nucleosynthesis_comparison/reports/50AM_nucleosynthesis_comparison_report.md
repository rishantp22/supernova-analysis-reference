# Comprehensive scientific and computational audit of 50AM nucleosynthesis

Audit date: 2026-10-18. Generated from completed measurements in `/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/`.

## Finding and publication assessment

**The differences are not merely formatting or row order. Neither plotting dataset currently supports an unqualified quantitative-publication claim.** Verified contributors include legacy absent-result zero filling, failure-associated legacy solver flags, different initial compositions, incorrectly labeled new trajectory time/density, mixed new output endpoints, and inconsistencies between selected last-time chunks and independently archived final abundances. The relative contribution of reaction-rate changes, NSE treatment and other nuclear physics cannot be isolated without controlled calculations using validated common inputs.

The all-particle statistics are complete for the available inputs. This does not mean every historical executable, input revision or multi-terabyte raw archive was recovered. Verified findings, evidence-supported mechanisms and unresolved run provenance are separated below. Original files were not modified, and no nuclear-network calculation was executed.

## Files, byte identity, rows and species

| Dataset | Full path | Bytes | Data rows | Valid unique IDs |
| --- | --- | ---: | ---: | ---: |
| Previous | `/data/pyoung3/supernova_data/50Am/analysis/50Am_plotting.out` | 748,723,480 | 945,628 | 945,628 |
| New | `/scratch/rpraka15/50Am_plotting_all_isotopes.out` | 1,766,352,794 | 937,228 | 937,227 |
| Intermediate export | `/scratch/rpraka15/final_chunks.out` | 1,438,981,268 | 937,229 | See prior intermediate-file audit |

Byte identity was ruled out by unequal file lengths before numerical comparison. SHA-256 is `7e88d0e59fbe27e9752837327ac3979a8ad5e1bec9c45aee4dab593e579436c6` for the previous file and `a5ad5ce7dc7a4beb8e75292fa0a4c64c56d71dde4bc16c8cedc223d554a55f24` for the new file. The alternate legacy copy under `/data/pyoung3/gsvance/supernova_data/50Am/analysis/50Am_plotting.out` has the same hash as the previous file. The requested literal name `50AmPlotting.out` was not found in the inventoried group paths.

The full ID join contains 937,227 valid shared IDs, 8,401 old-only IDs and zero new-only valid IDs. Both files have zero duplicate valid IDs, zero wrong-width rows and zero blank data lines. Ordering differs; old/new ID-inversion counts are 0/0. No comparison relies on row order. The new first data row contains invalid token `y193`, and ID 23 is absent relative to the complete chunks. Invalid records are preserved in audit tables, not silently corrected.

The original intermediate `final_chunks.out` is 1.44 GB decimal, not 19.7 GB. Its two known malformed physical tokens are `-0.045y2421` at ID 10016/vz and `7.02723e-0y5` at ID 10034/ax. The new exporter obtains physical values from the old plotting file, so those intermediate physical-token errors do not explain isotope differences. Prior streaming intermediate comparison results are retained at `/scratch/rpraka15/file_comparison_20261018/results-66028953/`; the current full chunk audit independently verifies all valid new isotope values.

Old/new files have 100/298 columns. There are 65/223 nuclear-species columns including neutron, 64 shared labels, but only 63 shared species with numeric new values because new neutron is entirely `NA`. New has 159 unique species labels; old uniquely exports 48Cr. Old/new elemental-total counts are 19/59; all new elemental totals are `NA`. Canonical names and original `X_{...}` labels are preserved in the complete inventory. Missing column, explicit NA and printed zero are never equivalent.

The raw legacy network stores 524 nuclides. The complete final EOC archive has 264 distinct (Z,A) rows and 245 species with any finite X, extending through 92Mo. The 13,396-entry current PRISM extent and the long heavy-isotope plotting header do not prove abundant heavy-species production. Many printed heavy-isotope zeros are selection/censoring outcomes. `species_inventory.csv`, `old_raw_species_inventory.csv`, `eoc_species_inventory.csv`, `initial_species_inventory.csv` and `prism_default_extent_species.csv` provide the complete inventories.

## Method and computational execution

Files were streamed into column-major float64 memory maps; IDs were joined explicitly. One species at a time supplies full-population differences, quantiles, extrema, relative changes, new/old ratios, positive-pair logarithmic differences, sign counts, newly zero/nonzero counts and yield contributions. No abundance floor or arbitrary logarithmic pseudocount is used. Exact equality is reported separately from approximate equality (rtol=1e-5, atol=0 for abundances; rtol=1e-9 for copied physical fields). Undefined ratios remain NaN. See the numerical validation report for rounding, ratio, trajectory and chunk-specific tolerances.

Integrated yields are `sum(mass_g * X) / 1.9889e33`, computed on documented cohorts. The shared plotting population represents 8.95708201662 solar masses in both files. No additional unbound-particle, fallback or decay selection is inferred. These are sums of available file populations, not certified physical ejected yields. The largest initial population process used about 3.23 GiB resident memory; other stages used memory maps, streamed archives and bounded batches. Heavy processing ran through Slurm; logs and accounting are retained.

The complete metadata inventory includes Sarah Vierling's scratch and accessible pyoung3 trees. Source content search covered 3,243 small candidate programs. The full compact final archive, all eight legacy HDF5 files, all 1,876 chunk files, all 943,900 candidate trajectories and all 943,900 seed files were processed. The ten multi-terabyte run archives were inspected only through bounded representative prefixes; their full run-by-run integrity remains unresolved.

## Verified causes and quantitative evidence

### Legacy absent-result zero filling

All 668,118 old plotting IDs without an available legacy HDF5 result have zero in every exported nuclear-species column. The exporter explicitly sums an empty query list to zero. Of shared plotting IDs, 668,114 such particles represent 6.46772831413 solar masses, whereas 269,113 with available raw legacy results represent 2.48935370249 solar masses. This produces apparent zero-to-positive changes independent of any nuclear rate change.

Legacy tooling separately computes aggregate yields with an unburned-material contribution. Its plotting-file path does not insert those per-particle unburned abundances. Therefore yields recomputed from this plotting table must not be equated with a separate historical corrected total-yield file. `yield_coverage_decomposition.csv` quantifies each species' contribution from the absent-raw cohort, and `old_raw_vs_new_selected_species.csv` compares the restricted common-raw cohort.

Additionally, 6,003,121 positive raw abundances are exported as zero among 17,760,640 compared cells for 64 legacy species. Omitted values reach 9.99986549343406e-7. Every other paired cell agrees within the query's seven-significant-digit precision (rtol=5.1e-7). This is verified species-dependent censoring, distinct from missing-particle zero filling. The old plotting 10Be column has no counterpart in the stored legacy network. The audit does not introduce its own truncation threshold.

### Legacy stored solver and history problems

Among 277,510 raw legacy records, 277,507 store looperror 1337 and three store 100. Available Fortran source assigns 1337 to failure branches and copies abundance/time results before checking the flag. This is strong evidence requiring failure-status resolution, not proof of the exact failing numerical branch or historical executable version. Raw closure reaches 1.38239028 at ID 157594. All 3,052,610 stored history times align with physical SDF times; 3,052,368 temperature/density pairs agree, with 242 inconsistent points across 22 particles. Examples include implausibly extreme and nonfinite stored values.

### New input trajectories and initial compositions

All 943,900 corrected-archive trajectories contain 61 points. Every temperature matches original SDF temperature/1e9 within printed precision. Every density matches native SDF density under the documented rounding allowance, despite its g/cm³ header. The conversion to cgs is 5.911824660348532e-6, corresponding to a factor about 169,152 between unconverted and cgs values. Every time label is a snapshot iteration number; the endpoint 2000 corresponds to physical time 157,558 s, with a nonlinear discrepancy throughout the series.

Ten archived runtime prefixes reproduce this input behavior, and each prefix's last saved abundances matches its new plotting row within printed precision. The inference is therefore anchored to actual calculations for those IDs. The full historical input revision for every tracer is not pinned by the available mutable control paths. Current source holds endpoint conditions constant beyond the candidate history endpoint, making a 166,000 s stop after a 2,000-labeled trajectory especially consequential; this is supported by source and inspected runtime histories, not inferred solely from file dates.

The complete initial-composition comparison found 943,900 of 943,900 particles differing in at least one shared native species at rtol=1e-5. The native hydro network is approximate (20 physical species), while PRISM seeds include a more detailed inventory, so small differences alone are not classified as errors. However, ID 917342 starts natively with X(1H)=0.66475999 and X(4He)=0.32527521, while its archived PRISM seed has X(16O)=0.945 before normalization. The corresponding archived runtime starts from that oxygen-rich composition at T≈0.000770219 GK. Such a seed change exists before burning and cannot be explained by a larger plotting header. Its producer/mapping must be recovered and validated. ID 191136 also differs materially in O/Si/S seed composition. Complete shared-native-species statistics and runtime examples are retained.

Reported progenitor Ye was independently checked: all old values match the native initial `Y_el` field within export precision. Thus the 495,014 zeros among valid new IDs originate in that native field, even though the exporter also has a generic missing-row zero fallback. Native composition-derived Ye disagrees with the stored field for much of the outer population. Neither copied zero Ye nor a partial isotope sum is a valid neutron-rich classification.

### Extraction and epoch inconsistencies

Every chunk grid and every valid new plotting isotope value passes the full comparison. There is no evidence of an ID or species-index shift in the final aggregation on valid rows. This locates remaining new-output differences upstream of final plotting assembly.

The independent EOC archive contains 938,732 unique IDs, 937,033 shared with valid new plotting IDs. There are 194 new valid IDs without EOC and 1,699 EOC IDs absent from valid new plotting IDs. Across 72,655,543 finite paired cells, 8,456 disagree at rtol=5.1e-6, affecting 176 particles. Another 5,682 positive plotting cells lack a finite EOC X. All finite X/Y pairs pass X=A Y at rtol=1e-12.

All 176 mismatching particles have chunk last times at or before 2000 s. ID 7753 has selected sum 1.977369060913544 versus EOC full sum 1.0000008545078554. Its 48 discrepant finite species are approximately doubled: 54Fe is 1.05158 in the plotting file versus 0.525792049882292 in EOC. The producer sums identical printed-time keys, which can merge repeated blocks. That mechanism, incomplete xtime output, restarts and stale final derivatives are supported possibilities; their individual causal assignments remain unresolved without the affected original 50AM time series and immutable run identities.

## All-isotope differences and yields

The following selected systems summarize the complete all-species tables; they do not limit the numerical scope. The ratios combine coverage, seed, physics, epoch and extraction effects and must not be described as pure improvements in nucleosynthesis.

| species | old_yield_shared_msun | new_yield_shared_msun | yield_ratio_shared | significant |
| --- | --- | --- | --- | --- |
| 24Mg | 0.003621552 | 0.1590208 | 43.90959 | 937227 |
| 25Mg | 6.144587e-05 | 0.001071889 | 17.44445 | 936006 |
| 26Mg | 0.005801344 | 0.004406169 | 0.7595083 | 934267 |
| 26Al | 1.251928e-07 | 3.704097e-06 | 29.58715 | 930661 |
| 27Al | 0.0005931054 | 0.004492307 | 7.574213 | 937225 |
| 40Ca | 0.006283636 | 0.008919129 | 1.419422 | 937227 |
| 41Ca | 1.040836e-06 | 2.419322e-05 | 23.24403 | 937222 |
| 42Ca | 2.721235e-05 | 0.001035561 | 38.05481 | 937227 |
| 43Ca | 1.367311e-07 | 5.098794e-06 | 37.29067 | 937227 |
| 44Ca | 5.066867e-08 | 9.935041e-06 | 196.0786 | 935115 |
| 44Ti | 7.578865e-05 | 7.039815e-06 | 0.09288745 | 937181 |
| 46Ti | 4.58938e-05 | 0.000170008 | 3.704378 | 937227 |
| 47Ti | 5.343591e-07 | 2.957218e-06 | 5.534141 | 937227 |
| 48Ti | 3.549205e-06 | 1.362208e-05 | 3.838066 | 937227 |
| 49Ti | 3.133361e-07 | 1.528538e-06 | 4.878269 | 937221 |
| 50Ti | 2.816634e-07 | 2.042952e-06 | 7.253168 | 927493 |
| 50Cr | 0.0001067298 | 0.0001366053 | 1.279916 | 937227 |
| 52Cr | 2.460978e-05 | 0.0002151125 | 8.740936 | 937226 |
| 53Cr | 2.340463e-06 | 3.617634e-06 | 1.545691 | 937219 |
| 54Cr | 7.377625e-07 | 1.520872e-05 | 20.61466 | 936412 |
| 54Fe | 0.001736217 | 0.04579191 | 26.37454 | 937227 |
| 56Fe | 0.00154804 | 0.005287405 | 3.415548 | 937215 |
| 57Fe | 5.303079e-05 | 0.0003224332 | 6.080113 | 937227 |
| 58Fe | 1.507545e-05 | 0.0006001493 | 39.80971 | 937218 |
| 60Fe | 1.400337e-07 | 6.502255e-08 | 0.464335 | 915631 |
| 56Ni | 0.05019201 | 0.05859997 | 1.167516 | 468100 |
| 58Ni | 0.001439428 | 0.02305692 | 16.01811 | 937227 |
| 60Ni | 7.174296e-05 | 0.0001932929 | 2.694243 | 937227 |
| 61Ni | 1.627909e-06 | 8.323263e-06 | 5.112854 | 937227 |
| 62Ni | 4.308534e-06 | 3.900747e-05 | 9.053537 | 937227 |
| 64Ni | 1.06676e-06 | 3.597278e-06 | 3.372151 | 916236 |


46Ca/48Ca are absent from the old plotting table; available raw-network comparisons on the 269,113-ID raw/new cohort are:

| species | finite_pairs | old_yield_shared_msun | new_yield_shared_msun | yield_ratio_shared | significant |
| --- | --- | --- | --- | --- | --- |
| 46Ca | 269113 | 1.105279e-08 | 1.991445e-08 | 1.801759 | 269113 |
| 48Ca | 269113 | 1.901743e-07 | 8.216139e-08 | 0.4320321 | 269113 |


The largest absolute particle-mass-weighted changes across all shared isotopes are:

| species | absolute_particle_yield_change_msun | max_abs_id | max_abs_old | max_abs_new |
| --- | --- | --- | --- | --- |
| 16O | 6.765141 | 35526 | 0.9521122 | 1.22969e-06 |
| 28Si | 0.5069348 | 32472 | 0.00322233 | 0.712299 |
| 4He | 0.258865 | 274944 | 0.38253 | 1.20966e-07 |
| 32S | 0.2375575 | 45519 | 0.000325802 | 0.421521 |
| 24Mg | 0.1576397 | 66397 | 0.001142204 | 0.0550343 |
| 12C | 0.09568121 | 210811 | 0.07450702 | 2.49879e-06 |
| 56Ni | 0.0912868 | 29053 | 1.901111e-05 | 0.861439 |
| 20Ne | 0.04499331 | 43589 | 0.0642779 | 5.50325e-08 |
| 54Fe | 0.04414279 | 7753 | 9.134148e-07 | 1.05158 |
| 29Si | 0.03837743 | 49491 | 4.506054e-05 | 0.00828876 |


The complete tables contain every shared isotope, all missingness states, quantiles, exact/approximate/significant counts, largest discrepancies and IDs. Separate all-available-population yields avoid silently dropping unmatched valid old particles. The extended raw/EOC audit compares all 262 shared raw species across 270,706 common raw/EOC IDs, retaining per-species finite availability.

## Physical interpretation and grain-analysis implications

All-species statistics are binned against reported peak temperature, density, Ye, mass and spatial radius, and against independently reconstructed physical cooling time. Figures use full-population two-dimensional histograms/hexbin counts where possible. Representative tracers were selected from population extremes, integrity failures or explicitly bounded archive prefixes, not used to estimate full-population statistics.

Those descriptive trends do not establish NSE, complete/incomplete silicon burning, alpha-rich or neutron-rich freeze-out, or oxygen/neon/carbon burning as the cause of a given discrepancy. Such regime assignments require valid time-dependent inputs, appropriate initial composition and trustworthy evolving Ye. The identified unit and seed problems change reaction exposure before subtle rate-library differences can be isolated.

Mg–Al ratios, radiogenic 26Mg interpretations, 48Ca–50Ti–54Cr correlations, Cr three-ratio diagrams and Fe-group predictions all require re-evaluation. The 26Al model abundance or number ratio is not an inferred initial Solar System ratio. No decay-corrected composition is assumed, and 26Al ground/isomer contributions cannot be separated from the available output keys.

The existing `/scratch/rpraka15/ccsn-presolar-oxide-analysis` code uses Cr isotope comparisons and Cr–Ti joint matching, with specific grain selection and lower-limit handling. These model-derived ratios, candidate tracer selections, likelihood/ranking results and publication figures need rerunning after validated replacement data exist. The configured observational workbook was not available in the inspected project inputs, so no grain-match result or observational agreement was fabricated. Existing code and research files were left unchanged.

## Answers to the requested questions

| question | answer |
| --- | --- |
| 1. Identical hydrodynamic trajectories? | Common SNSPH snapshot origin is numerically supported; the new candidate/runtime input time and density conventions differ. Copied plotting metadata alone was not used as proof. |
| 2. What changed? | Available workflows change Burn to PRISM, seed compositions, trajectory preparation, output epochs, species selection and censoring. Legacy missing-result zero filling and new extraction inconsistencies also contribute. |
| 3. Larger network or more exports? | Both representation changes are present: 524 stored legacy species versus a 13,396-species current PRISM extent, but only 222 numeric new plotting species. Effective historical active inventories are not established by column counts. |
| 4. Why different abundances despite equal tracer fields? | Those fields were copied. Independent input inspection finds altered seeds and thermodynamic histories, while old abundance coverage is incomplete. |
| 5. Which changed most? | 16O, 28Si, 4He and 32S dominate the absolute mass-weighted differences; complete all-isotope rankings and maximum-ID examples are provided. |
| 6. Mg–Al effects? | Shared-cohort 24Mg yield changes by 43.91×, 26Al by 29.59×, and 26Mg by 0.7595×. These combined differences cannot be treated as an isolated 26Al decay correction. |
| 7. Ca, especially 48Ca? | All available Ca species are audited. 48Ca and 46Ca are absent from the old plotting inventory but available in the raw network; their restricted-cohort comparisons are provided without substituting missing columns with zero. |
| 8. Ti, Cr, Fe, Ni and all others? | Every shared plotting isotope and every available raw-network species is included. Tables below and complete CSV/NPZ files quantify them; species unavailable in one output remain flagged. |
| 9. Burning regimes or Ye associations? | All-species population bins and physical cooling diagnostics exist. Input-unit, seed and reported-Ye problems prevent validated assignments to NSE/alpha-rich/neutron-rich burning from these comparisons alone. |
| 10. Physical consistency? | New full EOC closure is good; selected plotting closure fails for ID 7753. Legacy raw closure and history failures are documented. Stored hydro Ye is unsuitable for much of the population. |
| 11. Decay or output epochs? | Epochs differ and can contribute; no decay-only cause for the entire discrepancy is established. Incomplete chains and state information prevent complete corrections. |
| 12. Indexing, chunks or assembly? | Complete chunk grids and exact chunk-to-new isotope agreement rule out a final indexing/assembly mismatch on valid IDs. Upstream last-time/final-state inconsistency and approximately doubled ID 7753 abundances remain. |
| 13. Integrated yields? | All isotope yields, signed/absolute changes and discrepant-particle contributions are supplied on explicit cohorts. They are file-population sums, not certified ejecta yields. |
| 14. Scientific interpretations? | Ca–Ti–Cr and Mg–Al predictions change materially. The existing grain matching and correlated-isotope interpretations require recalculation after input/output validation; no observational fit was asserted. |
| 15. Which dataset is better validated? | Neither plotting file is publication-ready. The archived new EOC has stronger numerical closure, but this does not repair demonstrated input problems. The available old outputs have serious coverage and failure-status limitations. |
| 16. Further necessary tests? | Recover immutable run inputs/builds, validate seeds and physical histories, rebuild extraction from single identified final states, and perform a controlled small paired rerun before any larger production rerun. No network rerun was launched. |


## Deliverables and remaining tests

- [Key findings](KEY_FINDINGS.md).
- [Network and processing provenance](NUCLEAR_NETWORK_AND_PROCESSING_PROVENANCE.md).
- [Complete chunk audit](CHUNK_INTEGRITY_AUDIT.md).
- [Numerical validation](NUMERICAL_VALIDATION_REPORT.md).
- [Database schema and usage](../tables/README.md), including all matched-particle isotope records and stage-coverage flags.
- All-isotope yield/statistical/ratio tables in `../tables/`; 113 PNG and 113 vector PDF figures in `../figures/`.
- Reproduction programs and Slurm files in `../scripts/`, execution logs and resource accounting in `../logs/` and `../provenance/`.

Before publication: recover seed-generation provenance and per-run controls/builds; validate input time, density, initial composition and Ye against the intended hydro/stellar epochs; resolve failed legacy runs and restore explicit missing/unburned status; reconstruct new abundances from one identified final state per ID without repeated-time accumulation; reconcile all missing IDs and omitted species; choose a common physical output/decay epoch; then perform a controlled small paired calculation with fixed inputs and one configuration change at a time. A full production network rerun was not launched. These are remaining scientific requirements, not completed results.
