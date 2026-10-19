#!/usr/bin/env python3
"""Build reviewable scientific deliverables from completed audit measurements."""
from common import *
import datetime,os
def J(name):return json.loads((ROOT/'provenance'/f'{name}.json').read_text())
def C(name):return pd.read_csv(ROOT/'tables'/f'{name}.csv')
def md(df,cols=None):
    if cols is not None:df=df[cols]
    def cell(v):
        if pd.isna(v):return 'unavailable'
        if isinstance(v,(float,np.floating)):return f'{v:.7g}'
        return str(v).replace('|','/')
    return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(cell(x) for x in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
def write(name,text):
    (ROOT/name).write_text(text.strip()+'\n')
def ledger():
    seed=pd.read_csv(ROOT/'tables/initial_composition_population.csv.gz',usecols=['id']).id.to_numpy();traj=pd.read_csv(ROOT/'tables/trajectory_population.csv.gz',usecols=['id']).id.to_numpy()
    arrays={'old_plot':table('old')[2],'new_plot':table('new')[2],'chunks':np.load(ROOT/'tables/all_chunk_particles.npz')['id'],'legacy_raw':np.load(ROOT/'tables/raw_old_population_metadata.npz')['id'],'new_eoc':np.load(ROOT/'tables/eoc_particle_validation.npz')['id'],'archived_seed':seed,'archived_trajectory':traj}
    ids=np.unique(np.concatenate([v[v>=0] for v in arrays.values()]));d=pd.DataFrame({'id':ids});counts={}
    for k,v in arrays.items():d[k]=np.isin(ids,v[v>=0]);counts[k]=int(d[k].sum())
    d.to_csv(ROOT/'tables/tracer_stage_coverage.csv.gz',index=False,compression='gzip');save('provenance/stage_coverage_summary.json',dict(union_ids=len(ids),counts=counts,malformed_new_ID='y193 remains unassigned; no silent ID repair'))

def main():
    ledger();old=J('old_scan');new=J('new_scan');pop=J('population_summary');raw=J('raw_audit_summary');chunk=J('chunk_summary');eoc=J('eoc_summary');traj=J('trajectory_summary');cov=J('coverage_attribution');seed=J('seed_attribution_summary');phys=J('physical_source_audit');hist=J('legacy_history_vs_hydro');group=J('full_group_inventory_summary');source=J('full_group_source_leads');ye=J('ye_interpretation');alt=J('alternate_legacy_copy')
    stats=C('species_statistics');nuc=stats[stats.kind.eq('isotope')];inv=C('species_inventory');yields=C('yields_all_available_particles');einv=C('eoc_species_inventory');rawinv=C('old_raw_species_inventory');coverage=C('yield_coverage_decomposition');legacychecks=C('legacy_history_vs_hydro');bad=legacychecks[(legacychecks.matched_temperatures<legacychecks.points)|(legacychecks.matched_cgs_densities<legacychecks.points)];bad.to_csv(ROOT/'tables/legacy_history_inconsistent_particles.csv',index=False)
    interest=['24Mg','25Mg','26Mg','26Al','27Al','40Ca','41Ca','42Ca','43Ca','44Ca','44Ti','46Ti','47Ti','48Ti','49Ti','50Ti','50Cr','52Cr','53Cr','54Cr','53Mn','55Mn','54Fe','56Fe','57Fe','58Fe','60Fe','56Co','57Co','59Co','56Ni','57Ni','58Ni','60Ni','61Ni','62Ni','64Ni']
    focus=nuc[nuc.species.isin(interest)][['species','old_yield_shared_msun','new_yield_shared_msun','yield_ratio_shared','significant']]
    allfinite=nuc[nuc.finite_pairs>0];most=nuc.nlargest(10,'absolute_particle_yield_change_msun')[['species','absolute_particle_yield_change_msun','max_abs_id','max_abs_old','max_abs_new']]
    raw48=C('old_raw_vs_new_selected_species');raw48=raw48[raw48.species.isin(['46Ca','48Ca'])]
    figcount=len(list((ROOT/'figures').rglob('*.png')));pdfcount=len(list((ROOT/'figures').rglob('*.pdf')))
    provenance=f'''# Nuclear-network and processing provenance

## Evidence and search scope

The complete recursive metadata traversal covered {group['counts']['file']:,} files, {group['counts']['directory']:,} directories and {group['counts']['symlink']} symlinks under `/scratch/svierlin`, `/scratch/pyoung3`, `/data/pyoung3` and `/home/pyoung3`. `/home/pyoung3` was inaccessible; symlink targets were recorded without traversing them. A subsequent content search examined {source['source_files_examined']:,} small Python, shell, notebook and Fortran files. These are explicit scopes, not a claim to have inspected the contents of every group file. Full inventories, errors, source fingerprints and match locations are retained in `provenance/`.

The two legacy plotting copies under `/data/pyoung3/supernova_data/50Am/analysis/` and `/data/pyoung3/gsvance/supernova_data/50Am/analysis/` have identical SHA-256 hashes. Their inodes differ. A separate `50AmP` PRISM workflow also exists; its teammate scripts reference the unavailable historical `/home/gsvance/prismtools` path. It is not assumed to be Sarah's run. Direct ID 7753 files found in `jet3b_outputs` belong to another model and were not substituted for 50AM data.

The published model description identifies 50Am with 15S and describes TYCHO progenitors, neutrino-driven explosion treatment, SNSPH hydrodynamics and Burn postprocessing. This supports model nomenclature; it does not authenticate the particular files examined here. [Schulte et al., model methods](https://arxiv.org/pdf/2011.07459).

## Reconstructed processing chains

```mermaid
flowchart TD
 H[50AM SNSPH snapshots] --> B[Legacy histories in eight Burn HDF5 files]
 B --> Q[burn_query: selected and thresholded isotope records]
 Q --> O[dm_postprocess: old plotting file]
 H --> P[Candidate PRISM trajectories and initial compositions]
 P --> N[PRISM run archives]
 N --> X[xtime last-time aggregation: 222 selected isotopes]
 X --> C[938 A chunks and 938 weighted B chunks]
 C --> F[final_chunks.out]
 F --> A[New plotting abundances]
 O --> M[Physical fields joined by particle ID]
 M --> A
 N --> E[Independent final x/y extraction: EOC CSV archive]
```

The original hydro simulation and tracer extraction build are not fully pinned. The stored legacy histories and full candidate trajectory population were independently checked against the original SDF snapshots, providing numerical lineage beyond filenames or modification dates.

| Stage | Identified evidence | Qualification |
| --- | --- | --- |
| Native hydro | `/data/pyoung3/supernova_data/50Am/run3g_50Am6/`, 61 snapshots | Native initial/final physical fields checked over all 945,628 IDs |
| Legacy network | `/data/pyoung3/supernova_data/50Am/50Am.dir/out501.h5` through `out508.h5` | 277,510 unique stored IDs, identical 524-species Z/N arrays |
| Legacy source | `/data/pyoung3/burn/src/fnuclear.f`, `worker.f`; available `netsu.524`, `netweak.524`, `netwinv.524`, `nuc.in` and setup files | Current source/library availability does not prove the historical executable or input hashes |
| Legacy extraction | `/data/pyoung3/gsvance/data_mining/burn_query.c`, `run_query.py`, `dm_postprocess.py`, `particler.py`, `sn_utils.py`; original sbatch files | Query cutoff is configurable. An archived 26Al submission specifies `-a 12`; later query revisions and per-isotope cutoffs must not be inferred from that alone |
| New trajectories | `/scratch/svierlin/prism-v02/input/conditions/50Am_output.zip` and `50Am_output_tempcorr.zip` | Corrected archive comprehensively checked; controls point at a mutable unpacked path |
| New seeds | `/scratch/svierlin/prism-v02/input/conditions/50Am_abund.zip` | All 943,900 files inspected, closure/Ye and native shared species compared; original seed-generation script/run binding unresolved |
| PRISM code/config | `/scratch/svierlin/prism-v02/src/`, nine archived `control.json` files, current nuclear input files | Hashes retained; exact historical build and all-particle controls not recovered |
| PRISM outputs | Ten `eoc_batch_*.tar.gz` archives under `/data/pyoung3/groupdata/svierlin/50Am_data/` | Roughly 2.5 TB compressed; only bounded first-particle prefixes inspected, not full archive integrity |
| Independent final derivative | `/scratch/svierlin/prism-v02/EOC_csv.tar.gz`, produced by `EOC_compile.sh` | Entire archive processed: 938,732 IDs; producer skips existing CSVs, so stale/mixed revisions remain possible |
| Chunk producer | `/scratch/svierlin/prism-v02/sum_eoc_manifest_RP.sh` | Fixed 222-species selection; missing entries become printed zero; repeated time keys are summed |
| New exporter | `/scratch/rpraka15/final_chunks_export_20261018/make_plotting_format.py` | Physical fields copied from old file by ID; last-category chunk abundance used |

## Network and numerical settings

The legacy raw species inventory has 524 distinct (Z,A) entries. The currently available PRISM extent contains 13,396 permitted entries, Z=0 through 133. This is a configured extent, not a measurement of how many species were active in every historical calculation. All nine recovered controls share one nuclear-dataset configuration. Their 33 referenced current input files, including the extent, are readable and hashed. The active datasets reference Reaclib, Moller2003, NDI FRDM2012/AME2016 reaction and decay products, and NUBASE2016 decay/stability information. Paths and activation flags are preserved verbatim in `archived_control_records.json`.

Current PRISM `control.f90` defaults include initial timestep 1e-15 s, `minimum_mass_conservation=1e-8`, `nuc_mass_check=1e-3`, abundance threshold −1 and screening false unless overridden. These are source defaults, not a retrospectively verified solver configuration for every tracer. `progression.f90` normalizes initial X and computes evolving composition Ye for rate evaluation. `trajectory.f90` uses linear interpolation and holds endpoint temperature/density constant beyond the supplied interval. Archived controls specify stop times 2000 s for one inspected tracer and 166000 s for eight; the tenth prefix lacks a recovered control. Actual runtime headers corroborate the mixed endpoints.

The output schema identifies nuclides by Z/A, with no separate 26Al ground/isomer columns. No state-resolved 26Al ratio can be extracted from these tables. Exact legacy versus PRISM rate-library revisions, NSE prescriptions, partition functions, historical screening/weak/neutrino settings and solver builds remain insufficiently pinned for a controlled physics attribution. Existing source is evidence of available implementations, not proof that every feature was active. [PRISM-related methods paper](https://arxiv.org/abs/2008.06075).

The available PRISM output default lower bound is 1e-15, applied separately to X and Y. This explains why a CSV inventory row can have X or Y unavailable and why selected-species zero filling loses information. Original thresholds are part of output provenance, not a threshold applied by this audit. Where both are finite, all 111,953,173 EOC pairs satisfy X=A Y within relative tolerance 1e-12. In the legacy raw-to-plotting audit, all 17,760,640 paired cells for 64 species are explained by query rounding or positive-to-zero censoring. The latter affects 6,003,121 cells, with omitted values approaching 1e-6. The old 10Be plotting column has no counterpart in the stored legacy network; 10Be and 16N are present in the EOC inventory but absent from that legacy network.

## Earliest demonstrated divergence

1. The candidate PRISM seed compositions already differ from native initial SNSPH compositions before any new network integration. The native network is approximate, so modest species differences can be expected; large H/He-to-O differences require a seed-mapping explanation.
2. All corrected candidate histories preserve the hydro temperature sequence but substitute snapshot iteration for physical time and retain native density under a cgs label. Ten archived runtime prefixes numerically tie this behavior to actual new abundances.
3. Legacy raw records contain failure-associated status flags and {len(bad)} particles with inconsistent stored temperature/density histories. Missing raw results are then represented as zero in the old plotting output.
4. The selected new chunks reproduce the plotting file exactly, but 176 particles disagree with independent final EOC abundances. Their original 50AM time series and run identities are needed to distinguish incomplete histories, later restarts, stale final derivatives and repeated-time summation.

The same directory contains scripts for restarting from last saved abundances, replacing trajectory starts, and moving revised inputs. Their existence supports a revision risk but does not establish that a specific script caused a particular discrepancy. The inaccessible `eoc_manifest.txt` and absent unpacked `all_outputs/` prevent a full manifest-to-run binding. No group source or research data was altered.
'''
    write('reports/NUCLEAR_NETWORK_AND_PROCESSING_PROVENANCE.md',provenance)
    chunktext=f'''# Complete chunk integrity audit

All {chunk['chunks']} A/B chunk pairs were inspected, not sampled. A files total {chunk['A_bytes']:,} bytes and {chunk['A_rows']:,} records; B files total {chunk['B_bytes']:,} bytes and the same record count. The complete grid is 937,229 unique IDs × three categories × 222 species. Categories are nearest saved time to 100 s, nearest saved time to 1000 s, and last saved time; they are not guaranteed exact target epochs.

There are zero duplicate cross-chunk IDs, duplicate keys, missing ID/category/species cells, inconsistent A/B keys, negative/nonfinite abundances, or B mass-weighting discrepancies at relative tolerance 2e-5. A's column named `A` is Z and column `B` is nuclear mass number. B abundances use native particle-mass units of 1e-6 solar masses, not grams.

Every valid new plotting value exactly equals its last-category A abundance. The two chunk IDs absent from valid new-file IDs are 23 and 193; the latter has malformed token `y193` in the plotting file. No automatic repair was applied. 8,399 valid old IDs are absent from the chunks; see the stage ledger for differences against candidate inputs and independent EOC coverage.

Last saved times range from 0.936977 to 166000 seconds, with median 166000. The maximum last-category selected sum is 1.977369060913544 for ID 7753. There is no model/run identifier in a chunk row; complete grids alone cannot certify run identity or successful nuclear integration.

Independent archived final EOC comparisons found 8,456 finite mismatches in 176 particles using rtol=5.1e-6 and atol=0. All these particles have chunk endpoints at or before 2000 s. Particle 7753 has approximately doubled values for 48 finite discrepant species. A/B internal consistency is therefore compatible with an upstream extraction error. The producer's printed-time summation is a demonstrated mechanism capable of merging repeated blocks, but the actual affected 50AM time series were not recovered, so this mechanism remains an evidence-supported hypothesis.

During the first pass, 24 chunks changed filesystem metadata while being read. Their failed attempts are preserved. The retry required matching content hashes when metadata changed; all 938 pairs completed successfully. End-of-audit source stability checks rehashed changed-metadata sources and rehashed both plotting files. This is not evidence that the initial metadata changes altered abundance bytes. Receipts include paths, sizes, metadata and SHA-256.

Detailed machine-readable evidence: `tables/chunk_integrity.csv`, `tables/all_chunk_particles.npz`, `tables/eoc_mismatches_with_chunk_epoch.csv`, `tables/eoc_plot_mismatch_cells.csv`, `provenance/chunks/`, `provenance/source_stability.json`.
'''
    write('reports/CHUNK_INTEGRITY_AUDIT.md',chunktext)
    validation=f'''# Numerical and physical validation

## Population and physical quantities

The plotting comparison includes all 937,227 shared unique valid IDs. Every one of the 15 shared physical columns is exactly equal after numeric parsing, with no merely approximate or significant physical mismatch. This equality is constructed by the new exporter. An independent source check also confirms all 945,628 old-file IDs' final x/y/z, velocity, acceleration, mass, h and density against native SDF fields after conversion, at rtol=1e-5 and atol=0. Initial reported Ye also matches the native stored `Y_el` field at this precision.

Units: positions and h in cm; velocities cm/s; accelerations cm/s²; mass g; density g/cm³; peak temperature K; Ye dimensionless; isotope/element X dimensionless mass fractions. Native length, time and mass factors are 6.955e10 cm, 100 s and 1.9889e27 g. Legacy 'peak density' is density at the sampled peak-temperature point, not an independently maximized density. The old early-snapshot peak window differs from a maximum over all 61 snapshots.

**Correction to an early interpretation:** although the legacy exporter has a missing-record-to-zero path, the observed copied zero Ye values match zeros already in the original hydro `Y_el` field. Thus their presence is not attributed to an exporter lookup failure. There are {ye['reported_Ye_zero']:,} such zeros among valid new IDs. Native composition-derived Ye exposes their inconsistency; neither these zeros nor partial plotting Ye should be used to classify neutron-rich freeze-out.

All 3,052,610 legacy history times match physical snapshot times at rtol=1e-6, atol=1e-9 s. Temperature and converted density both match at 3,052,368 points with rtol=1e-10, atol=0. The remaining 242 points belong to {len(bad)} particles; their stored values include extreme signed magnitudes and one nonfinite point. Full particle flags and examples are retained, without repairing the HDF5 data.

## Mass closure and electron fraction

Available-isotope plotting sums exclude elemental totals and count each nuclide once. The inventories are incomplete, so deficits are not automatically errors. Over-unity sums are inconsistent with valid mass fractions even for incomplete inventories.

{md(C('closure_statistics'),['dataset','particles','isotope_columns','rows_missing','over_unity_1e5','sum_min','sum_p50','sum_max','max_sum_id'])}

All 524 legacy raw species give sum-X range {raw['raw_sumX']['min']:.10g}–{raw['raw_sumX']['max']:.10g}; {raw['raw_over_unity_1e5']} particles exceed 1.00001 and {raw['raw_under_unity_1e5']} fall below 0.99999. Maximum is ID 157594. These are checks of saved results, not a successful-solver certificate.

The full independent EOC archive gives finite-X sums from {eoc['sumX']['min']:.12g} to {eoc['sumX']['max']:.12g}; all 938,732 are within 1e-5 of unity. Its composition-derived Ye ranges {eoc['Ye']['min']:.10g}–{eoc['Ye']['max']:.10g}. Censored/missing X remains recorded. Closure cannot validate trajectory units, seeds, reaction rates or historical run identity.

Native initial composition uses slots 1–20; the source explicitly excludes duplicate free-nucleon slots 21–22. Native composition closes within approximately 8e-6, while native stored Ye disagrees severely in much of the outer population. Archived PRISM seeds have sums 0.999232410975–1.000814473939 and are normalized explicitly by the network initialization. The audit preserves original seed values and separately labels normalized diagnostics.

All new elemental totals (59 columns) are unavailable. Legacy elemental totals are compared with available isotope subtotals in `elemental_total_checks.csv`; a subtotal deficit can reflect unexported isotopes or thresholded queries and is not sufficient evidence of incorrect elemental mass. No new total is synthesized and substituted into the original file.

## Ratios, decay and censoring

Number ratios use `(X_i/A_i)/(X_j/A_j)`. The 172 audited pairs include every available isotope against a declared within-element reference plus all requested Mg–Al, Ca, Ti, Cr and Fe pairs. Missing columns and zero denominators remain explicit; no pseudocount or abundance floor is used. The all-particle ratio files and quantiles are in `tables/ratios/` and `number_ratio_statistics.csv`.

The old/new tables are not demonstrated to be at one common freeze-out or decay epoch. Chain subtotal checks test whether *isolated* decay alone could explain a difference; they do not apply decay corrections and do not prove decay when a sum agrees. The 44Ti→44Sc→44Ca chain is unavailable as a complete plotting chain because 44Sc is omitted. The 26Al output has no state separation. No inferred initial Solar System 26Al/27Al is equated with a model mass fraction or model number ratio.

{md(C('decay_chain_diagnostics'),['chain','old_missing','new_missing','valid_closed_chain_pairs','inconsistent_with_isolated_decay_1e5'])}

Thresholded source output is different from exact zero. The EOC comparison records 135,360,101 plotting-zero cells without a finite final-X counterpart and 5,682 plotting-positive cells without one. These cases cannot be forced into ordinary relative/log comparisons. Source X=A Y consistency holds for every finite pair tested.

## Tolerances and numerical limits

The primary abundance comparison uses rtol=1e-5, atol=0; copied physical fields rtol=1e-9, atol=0; isotope number ratios rtol=2e-5, atol=0. Exact matches are separately counted. Raw query rounding checks use 5.1e-7; EOC/chunk output rounding 5.1e-6; B weighting 2e-5. All are computational tolerances, not physical error bars. The candidate density export uses rtol=5.1e-6 plus atol=5.1e-6 native units to accommodate five printed decimal places; without that documented rounding allowance, small native densities appear discrepant.

All statistics use the full relevant population. Distribution percentiles omit undefined/nonfinite values and accompany availability/overflow counts. Log differences are evaluated only for positive pairs without clipping tiny abundances. No input normalization is imposed except the explicitly labeled diagnostic reproducing the network's own initial normalization.
'''
    write('reports/NUMERICAL_VALIDATION_REPORT.md',validation)
    questions=[
    ('Identical hydrodynamic trajectories?', 'Common SNSPH snapshot origin is numerically supported; the new candidate/runtime input time and density conventions differ. Copied plotting metadata alone was not used as proof.'),
    ('What changed?', 'Available workflows change Burn to PRISM, seed compositions, trajectory preparation, output epochs, species selection and censoring. Legacy missing-result zero filling and new extraction inconsistencies also contribute.'),
    ('Larger network or more exports?', 'Both representation changes are present: 524 stored legacy species versus a 13,396-species current PRISM extent, but only 222 numeric new plotting species. Effective historical active inventories are not established by column counts.'),
    ('Why different abundances despite equal tracer fields?', 'Those fields were copied. Independent input inspection finds altered seeds and thermodynamic histories, while old abundance coverage is incomplete.'),
    ('Which changed most?', '16O, 28Si, 4He and 32S dominate the absolute mass-weighted differences; complete all-isotope rankings and maximum-ID examples are provided.'),
    ('Mg–Al effects?', 'Shared-cohort 24Mg yield changes by 43.91×, 26Al by 29.59×, and 26Mg by 0.7595×. These combined differences cannot be treated as an isolated 26Al decay correction.'),
    ('Ca, especially 48Ca?', 'All available Ca species are audited. 48Ca and 46Ca are absent from the old plotting inventory but available in the raw network; their restricted-cohort comparisons are provided without substituting missing columns with zero.'),
    ('Ti, Cr, Fe, Ni and all others?', 'Every shared plotting isotope and every available raw-network species is included. Tables below and complete CSV/NPZ files quantify them; species unavailable in one output remain flagged.'),
    ('Burning regimes or Ye associations?', 'All-species population bins and physical cooling diagnostics exist. Input-unit, seed and reported-Ye problems prevent validated assignments to NSE/alpha-rich/neutron-rich burning from these comparisons alone.'),
    ('Physical consistency?', 'New full EOC closure is good; selected plotting closure fails for ID 7753. Legacy raw closure and history failures are documented. Stored hydro Ye is unsuitable for much of the population.'),
    ('Decay or output epochs?', 'Epochs differ and can contribute; no decay-only cause for the entire discrepancy is established. Incomplete chains and state information prevent complete corrections.'),
    ('Indexing, chunks or assembly?', 'Complete chunk grids and exact chunk-to-new isotope agreement rule out a final indexing/assembly mismatch on valid IDs. Upstream last-time/final-state inconsistency and approximately doubled ID 7753 abundances remain.'),
    ('Integrated yields?', 'All isotope yields, signed/absolute changes and discrepant-particle contributions are supplied on explicit cohorts. They are file-population sums, not certified ejecta yields.'),
    ('Scientific interpretations?', 'Ca–Ti–Cr and Mg–Al predictions change materially. The existing grain matching and correlated-isotope interpretations require recalculation after input/output validation; no observational fit was asserted.'),
    ('Which dataset is better validated?', 'Neither plotting file is publication-ready. The archived new EOC has stronger numerical closure, but this does not repair demonstrated input problems. The available old outputs have serious coverage and failure-status limitations.'),
    ('Further necessary tests?', 'Recover immutable run inputs/builds, validate seeds and physical histories, rebuild extraction from single identified final states, and perform a controlled small paired rerun before any larger production rerun. No network rerun was launched.')]
    report=f'''# Comprehensive scientific and computational audit of 50AM nucleosynthesis

Audit date: 2026-10-18. Generated from completed measurements in `/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/`.

## Finding and publication assessment

**The differences are not merely formatting or row order. Neither plotting dataset currently supports an unqualified quantitative-publication claim.** Verified contributors include legacy absent-result zero filling, failure-associated legacy solver flags, different initial compositions, incorrectly labeled new trajectory time/density, mixed new output endpoints, and inconsistencies between selected last-time chunks and independently archived final abundances. The relative contribution of reaction-rate changes, NSE treatment and other nuclear physics cannot be isolated without controlled calculations using validated common inputs.

The all-particle statistics are complete for the available inputs. This does not mean every historical executable, input revision or multi-terabyte raw archive was recovered. Verified findings, evidence-supported mechanisms and unresolved run provenance are separated below. Original files were not modified, and no nuclear-network calculation was executed.

## Files, byte identity, rows and species

| Dataset | Full path | Bytes | Data rows | Valid unique IDs |
| --- | --- | ---: | ---: | ---: |
| Previous | `{old['source']['path']}` | {old['source']['bytes']:,} | {old['rows']:,} | {old['unique_valid_ids']:,} |
| New | `{new['source']['path']}` | {new['source']['bytes']:,} | {new['rows']:,} | {new['unique_valid_ids']:,} |
| Intermediate export | `/scratch/rpraka15/final_chunks.out` | 1,438,981,268 | 937,229 | See prior intermediate-file audit |

Byte identity was ruled out by unequal file lengths before numerical comparison. SHA-256 is `{old['sha256']}` for the previous file and `{new['sha256']}` for the new file. The alternate legacy copy under `/data/pyoung3/gsvance/supernova_data/50Am/analysis/50Am_plotting.out` has the same hash as the previous file. The requested literal name `50AmPlotting.out` was not found in the inventoried group paths.

The full ID join contains {pop['shared_unique_valid_ids']:,} valid shared IDs, {pop['old_only_ids']:,} old-only IDs and zero new-only valid IDs. Both files have zero duplicate valid IDs, zero wrong-width rows and zero blank data lines. Ordering differs; old/new ID-inversion counts are {old['id_inversions']:,}/{new['id_inversions']:,}. No comparison relies on row order. The new first data row contains invalid token `y193`, and ID 23 is absent relative to the complete chunks. Invalid records are preserved in audit tables, not silently corrected.

The original intermediate `final_chunks.out` is 1.44 GB decimal, not 19.7 GB. Its two known malformed physical tokens are `-0.045y2421` at ID 10016/vz and `7.02723e-0y5` at ID 10034/ax. The new exporter obtains physical values from the old plotting file, so those intermediate physical-token errors do not explain isotope differences. Prior streaming intermediate comparison results are retained at `/scratch/rpraka15/file_comparison_20261018/results-66028953/`; the current full chunk audit independently verifies all valid new isotope values.

Old/new files have {len(old['columns'])}/{len(new['columns'])} columns. There are 65/223 nuclear-species columns including neutron, 64 shared labels, but only 63 shared species with numeric new values because new neutron is entirely `NA`. New has 159 unique species labels; old uniquely exports 48Cr. Old/new elemental-total counts are 19/59; all new elemental totals are `NA`. Canonical names and original `X_{{...}}` labels are preserved in the complete inventory. Missing column, explicit NA and printed zero are never equivalent.

The raw legacy network stores 524 nuclides. The complete final EOC archive has {len(einv)} distinct (Z,A) rows and {(einv.finite_x>0).sum()} species with any finite X, extending through 92Mo. The 13,396-entry current PRISM extent and the long heavy-isotope plotting header do not prove abundant heavy-species production. Many printed heavy-isotope zeros are selection/censoring outcomes. `species_inventory.csv`, `old_raw_species_inventory.csv`, `eoc_species_inventory.csv`, `initial_species_inventory.csv` and `prism_default_extent_species.csv` provide the complete inventories.

## Method and computational execution

Files were streamed into column-major float64 memory maps; IDs were joined explicitly. One species at a time supplies full-population differences, quantiles, extrema, relative changes, new/old ratios, positive-pair logarithmic differences, sign counts, newly zero/nonzero counts and yield contributions. No abundance floor or arbitrary logarithmic pseudocount is used. Exact equality is reported separately from approximate equality (rtol=1e-5, atol=0 for abundances; rtol=1e-9 for copied physical fields). Undefined ratios remain NaN. See the numerical validation report for rounding, ratio, trajectory and chunk-specific tolerances.

Integrated yields are `sum(mass_g * X) / 1.9889e33`, computed on documented cohorts. The shared plotting population represents {pop['shared_particle_mass_old_msun']:.11f} solar masses in both files. No additional unbound-particle, fallback or decay selection is inferred. These are sums of available file populations, not certified physical ejected yields. The largest initial population process used about 3.23 GiB resident memory; other stages used memory maps, streamed archives and bounded batches. Heavy processing ran through Slurm; logs and accounting are retained.

The complete metadata inventory includes Sarah Vierling's scratch and accessible pyoung3 trees. Source content search covered {source['source_files_examined']:,} small candidate programs. The full compact final archive, all eight legacy HDF5 files, all 1,876 chunk files, all 943,900 candidate trajectories and all 943,900 seed files were processed. The ten multi-terabyte run archives were inspected only through bounded representative prefixes; their full run-by-run integrity remains unresolved.

## Verified causes and quantitative evidence

### Legacy absent-result zero filling

All 668,118 old plotting IDs without an available legacy HDF5 result have zero in every exported nuclear-species column. The exporter explicitly sums an empty query list to zero. Of shared plotting IDs, 668,114 such particles represent 6.46772831413 solar masses, whereas 269,113 with available raw legacy results represent 2.48935370249 solar masses. This produces apparent zero-to-positive changes independent of any nuclear rate change.

Legacy tooling separately computes aggregate yields with an unburned-material contribution. Its plotting-file path does not insert those per-particle unburned abundances. Therefore yields recomputed from this plotting table must not be equated with a separate historical corrected total-yield file. `yield_coverage_decomposition.csv` quantifies each species' contribution from the absent-raw cohort, and `old_raw_vs_new_selected_species.csv` compares the restricted common-raw cohort.

Additionally, 6,003,121 positive raw abundances are exported as zero among 17,760,640 compared cells for 64 legacy species. Omitted values reach 9.99986549343406e-7. Every other paired cell agrees within the query's seven-significant-digit precision (rtol=5.1e-7). This is verified species-dependent censoring, distinct from missing-particle zero filling. The old plotting 10Be column has no counterpart in the stored legacy network. The audit does not introduce its own truncation threshold.

### Legacy stored solver and history problems

Among 277,510 raw legacy records, 277,507 store looperror 1337 and three store 100. Available Fortran source assigns 1337 to failure branches and copies abundance/time results before checking the flag. This is strong evidence requiring failure-status resolution, not proof of the exact failing numerical branch or historical executable version. Raw closure reaches 1.38239028 at ID 157594. All 3,052,610 stored history times align with physical SDF times; 3,052,368 temperature/density pairs agree, with 242 inconsistent points across {len(bad)} particles. Examples include implausibly extreme and nonfinite stored values.

### New input trajectories and initial compositions

All 943,900 corrected-archive trajectories contain 61 points. Every temperature matches original SDF temperature/1e9 within printed precision. Every density matches native SDF density under the documented rounding allowance, despite its g/cm³ header. The conversion to cgs is 5.911824660348532e-6, corresponding to a factor about 169,152 between unconverted and cgs values. Every time label is a snapshot iteration number; the endpoint 2000 corresponds to physical time 157,558 s, with a nonlinear discrepancy throughout the series.

Ten archived runtime prefixes reproduce this input behavior, and each prefix's last saved abundances matches its new plotting row within printed precision. The inference is therefore anchored to actual calculations for those IDs. The full historical input revision for every tracer is not pinned by the available mutable control paths. Current source holds endpoint conditions constant beyond the candidate history endpoint, making a 166,000 s stop after a 2,000-labeled trajectory especially consequential; this is supported by source and inspected runtime histories, not inferred solely from file dates.

The complete initial-composition comparison found {seed['particles_any_shared_isotope_difference']:,} of {seed['ids']:,} particles differing in at least one shared native species at rtol=1e-5. The native hydro network is approximate (20 physical species), while PRISM seeds include a more detailed inventory, so small differences alone are not classified as errors. However, ID 917342 starts natively with X(1H)=0.66475999 and X(4He)=0.32527521, while its archived PRISM seed has X(16O)=0.945 before normalization. The corresponding archived runtime starts from that oxygen-rich composition at T≈0.000770219 GK. Such a seed change exists before burning and cannot be explained by a larger plotting header. Its producer/mapping must be recovered and validated. ID 191136 also differs materially in O/Si/S seed composition. Complete shared-native-species statistics and runtime examples are retained.

Reported progenitor Ye was independently checked: all old values match the native initial `Y_el` field within export precision. Thus the {ye['reported_Ye_zero']:,} zeros among valid new IDs originate in that native field, even though the exporter also has a generic missing-row zero fallback. Native composition-derived Ye disagrees with the stored field for much of the outer population. Neither copied zero Ye nor a partial isotope sum is a valid neutron-rich classification.

### Extraction and epoch inconsistencies

Every chunk grid and every valid new plotting isotope value passes the full comparison. There is no evidence of an ID or species-index shift in the final aggregation on valid rows. This locates remaining new-output differences upstream of final plotting assembly.

The independent EOC archive contains 938,732 unique IDs, 937,033 shared with valid new plotting IDs. There are 194 new valid IDs without EOC and 1,699 EOC IDs absent from valid new plotting IDs. Across 72,655,543 finite paired cells, 8,456 disagree at rtol=5.1e-6, affecting 176 particles. Another 5,682 positive plotting cells lack a finite EOC X. All finite X/Y pairs pass X=A Y at rtol=1e-12.

All 176 mismatching particles have chunk last times at or before 2000 s. ID 7753 has selected sum 1.977369060913544 versus EOC full sum 1.0000008545078554. Its 48 discrepant finite species are approximately doubled: 54Fe is 1.05158 in the plotting file versus 0.525792049882292 in EOC. The producer sums identical printed-time keys, which can merge repeated blocks. That mechanism, incomplete xtime output, restarts and stale final derivatives are supported possibilities; their individual causal assignments remain unresolved without the affected original 50AM time series and immutable run identities.

## All-isotope differences and yields

The following selected systems summarize the complete all-species tables; they do not limit the numerical scope. The ratios combine coverage, seed, physics, epoch and extraction effects and must not be described as pure improvements in nucleosynthesis.

{md(focus)}

46Ca/48Ca are absent from the old plotting table; available raw-network comparisons on the 269,113-ID raw/new cohort are:

{md(raw48,['species','finite_pairs','old_yield_shared_msun','new_yield_shared_msun','yield_ratio_shared','significant'])}

The largest absolute particle-mass-weighted changes across all shared isotopes are:

{md(most)}

The complete tables contain every shared isotope, all missingness states, quantiles, exact/approximate/significant counts, largest discrepancies and IDs. Separate all-available-population yields avoid silently dropping unmatched valid old particles. The extended raw/EOC audit compares all 262 shared raw species across 270,706 common raw/EOC IDs, retaining per-species finite availability.

## Physical interpretation and grain-analysis implications

All-species statistics are binned against reported peak temperature, density, Ye, mass and spatial radius, and against independently reconstructed physical cooling time. Figures use full-population two-dimensional histograms/hexbin counts where possible. Representative tracers were selected from population extremes, integrity failures or explicitly bounded archive prefixes, not used to estimate full-population statistics.

Those descriptive trends do not establish NSE, complete/incomplete silicon burning, alpha-rich or neutron-rich freeze-out, or oxygen/neon/carbon burning as the cause of a given discrepancy. Such regime assignments require valid time-dependent inputs, appropriate initial composition and trustworthy evolving Ye. The identified unit and seed problems change reaction exposure before subtle rate-library differences can be isolated.

Mg–Al ratios, radiogenic 26Mg interpretations, 48Ca–50Ti–54Cr correlations, Cr three-ratio diagrams and Fe-group predictions all require re-evaluation. The 26Al model abundance or number ratio is not an inferred initial Solar System ratio. No decay-corrected composition is assumed, and 26Al ground/isomer contributions cannot be separated from the available output keys.

The existing `/scratch/rpraka15/ccsn-presolar-oxide-analysis` code uses Cr isotope comparisons and Cr–Ti joint matching, with specific grain selection and lower-limit handling. These model-derived ratios, candidate tracer selections, likelihood/ranking results and publication figures need rerunning after validated replacement data exist. The configured observational workbook was not available in the inspected project inputs, so no grain-match result or observational agreement was fabricated. Existing code and research files were left unchanged.

## Answers to the requested questions

{md(pd.DataFrame([dict(question=f'{i}. {q}',answer=a) for i,(q,a) in enumerate(questions,1)]))}

## Deliverables and remaining tests

- [Key findings](KEY_FINDINGS.md).
- [Network and processing provenance](NUCLEAR_NETWORK_AND_PROCESSING_PROVENANCE.md).
- [Complete chunk audit](CHUNK_INTEGRITY_AUDIT.md).
- [Numerical validation](NUMERICAL_VALIDATION_REPORT.md).
- [Database schema and usage](../tables/README.md), including all matched-particle isotope records and stage-coverage flags.
- All-isotope yield/statistical/ratio tables in `../tables/`; {figcount} PNG and {pdfcount} vector PDF figures in `../figures/`.
- Reproduction programs and Slurm files in `../scripts/`, execution logs and resource accounting in `../logs/` and `../provenance/`.

Before publication: recover seed-generation provenance and per-run controls/builds; validate input time, density, initial composition and Ye against the intended hydro/stellar epochs; resolve failed legacy runs and restore explicit missing/unburned status; reconstruct new abundances from one identified final state per ID without repeated-time accumulation; reconcile all missing IDs and omitted species; choose a common physical output/decay epoch; then perform a controlled small paired calculation with fixed inputs and one configuration change at a time. A full production network rerun was not launched. These are remaining scientific requirements, not completed results.
'''
    write('reports/50AM_nucleosynthesis_comparison_report.md',report)
    # Full inventory unions make network omissions explicit, without zero filling.
    erows={int(r.Z)*512+int(r.A):r for r in einv.itertuples(index=False)};rset={int(r.Z)*512+int(r.A) for r in rawinv.itertuples(index=False)};extent=C('prism_default_extent_species');pset={int(r.Z)*512+int(r.A) for r in extent.itertuples(index=False)};plot_old={species(c)[1]*512+species(c)[2] for c in old['columns'] if species(c)[0]=='isotope'};plot_new={species(c)[1]*512+species(c)[2] for c in new['columns'] if species(c)[0]=='isotope'}
    allcodes=sorted(rset|pset|set(erows));csvsave('tables/complete_network_export_inventory.csv',[dict(Z=c//512,A=c%512,in_legacy_raw=c in rset,in_current_prism_extent=c in pset,in_final_eoc=c in erows,eoc_has_finite_X=bool(c in erows and erows[c].finite_x>0),in_old_plot=c in plot_old,in_new_plot=c in plot_new,new_plot_numeric=bool(c in plot_new and c!=1)) for c in allcodes])
    manifest=[]
    for base in ['scripts','reports']:
        for p in sorted((ROOT/base).glob('*')):
            if p.is_file():manifest.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)))
    csvsave('provenance/audit_code_report_manifest.csv',manifest);log('Scientific report, provenance, validation, chunk report and database ledger written')
if __name__=='__main__':main()
