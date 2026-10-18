# How to read the data and this audit

## The hierarchy

Think of each model as an experiment. A hydrodynamic **snapshot** is a photograph
of its fluid particles at one simulation time. An SPH **particle** is a moving
fluid parcel identified by `ident`. Its nucleosynthesis **trajectory** describes
temperature and density over time. PRISM follows nuclear reactions along that
trajectory and writes one output directory per particle.

```text
model
├── many SDF snapshots: one record per surviving particle
│   └── ident, position, velocity, acceleration, mass, h, rho, temp, ...
└── many PRISM particle directories: one reaction-network calculation each
    ├── control.json → trajectory and initial-composition paths
    ├── x.txt → rows (Z, A, final isotope mass fraction)
    ├── xz.txt → rows (Z, final element mass fraction)
    └── xtime.txt → time blocks, then rows (Z, A, isotope mass fraction)
```

The 50Am PRISM directories are packaged into batches of tar.gz archives. A batch
number is a storage partition, not a model, snapshot time or isotope. Archive
members such as `191136/x.txt` still identify the particle.

## Your requested fields

| Request | Source and interpretation |
|---|---|
| Id | SDF `ident`, numeric PRISM directory/member name; join within model |
| x, y, z | Snapshot position; convert code length using that run's factor |
| vx, vy, vz | Snapshot velocity; lengthCF / timeCF |
| ax, ay, az | Snapshot acceleration; lengthCF / timeCF² |
| mass | SPH parcel mass; multiply code mass by massCF for grams |
| h | SPH smoothing length, not a dust-grain radius |
| Density | Snapshot `rho`; massCF / lengthCF³ |
| Peak temp | Maximum over a validated particle history, not snapshot temperature |
| Peak density | Maximum over a validated particle density history |
| Y_{e} | Investigate stored `Y_el`; independently calculate sum(Z/A × X) for a chosen composition and epoch |
| X_{n} | Free-neutron mass fraction: PRISM Z=0, A=1 |
| 12C, 26Al, etc. | `x.txt` or selected `xtime.txt` block, match exact Z,A |
| C, Al, etc. | `xz.txt` elemental total over the whole network |

`Y` in the isotope/element list means **yttrium**, while `Y_{e}` is electron
fraction and PRISM's `y.txt` is abundance X/A. `xa.txt` groups by mass number and
must not be used to retrieve a unique isotope: nuclei with different Z can have
the same A. Do not double-count elemental totals and their constituent isotopes.

Temperature in PRISM time headers is GK; multiply by 10^9 for K. A filename
suffix is an iteration, not seconds. A maximum over saved history is only a
sampled maximum over that history's coverage, not a guaranteed continuous peak.
Distinguish full hydro history from PRISM's integration interval.

The complete machine-readable mapping is in scratch
`tables/requested_field_map.csv`, one row per model/request. `stored_samples`
counts inspected particles only. `requested_values_sample.csv` gives values,
statuses, selected end times and exact sources for those examples. Blank values
are intentional where a source is missing, a sparse row is omitted, or a
hydrodynamic epoch has not been chosen. The separate ID join demonstration
shows hydro values without pretending they are at the PRISM abundance epoch.

## What the plots show

1. **Snapshot times:** why filename suffixes cannot stand in for physical time.
2. **Particle geometry:** a systematic 20,000-record sample per chosen snapshot,
   colored by temperature. Panels use different end times and color scales;
   use them to inspect structure, not compare models at equal age.
3. **Requested composition coverage:** stored versus omitted requested cells.
   This is a bounded, nonrandom audit sample, not a population estimate.
4. **PRISM sanity:** normalization and actual saved endpoint differences.
5. **Input trajectories:** the inputs as written, with disputed units clearly
   labeled; a blank 50Am panel means the configured source path is missing.
6. **Hydro composition flags:** nonunit sums and stored zero electron fractions.
7. **Unit crosscheck:** empirical comparisons to the corresponding raw snapshots.

## Findings to resolve before producing a science table

The source-code and output checks in the linked provenance report are material:
time histories determine how long reactions act, and density affects their
rates. A corrected postprocessing run may be necessary if the historical run
consumed the disputed inputs. Do not change a final CSV's units and assume its
abundances are repaired. The audit preserves the original values and records
the discrepancy rather than editing master data.

The selected jet3b and vconv late snapshot samples store zero `Y_el`. A field's
presence in the schema does not prove it was populated meaningfully. The
sampled 50Am/jet3b hydro nuclear-slot sums reach approximately 1.03; diagnose
their network semantics and evolution before treating them as a normalized
composition. PRISM normalization looks much closer to unity in the inspected
examples, which does not validate its input physics.

The `vconv` PRISM folder, `vconv` binaries and `vconvL` text products have related
provenance, but their complete equivalence is unverified. Keep both the requested
label and the observed label in future exports. Initial-composition paths named
50Am are shared by other models and do not establish identical later evolution.

Old per-isotope CSV chunks select several saved epochs, insert some zeros, and
later drop time or particle identifiers. They cannot reconstruct the complete
requested per-particle dataset. The existing files are preserved for comparison.

## Contract for a full extraction

The next production table needs a key `(observed_model, particle_id, epoch)`.
For each value preserve source path, units, conversion-factor provenance,
quantity definition, and a status such as stored, derived, outside network,
missing source, sparse omission, or failed validation. Store particle mass
separately from isotope mass fraction; isotope mass is their product only when
model/particle/epoch and mass-unit assumptions are justified.

For an isotope absent from sparse output, do not manufacture zero. For elemental
mass fractions use xz, or a qualified sum over the full network at a requested
time. Choose whether radioactive abundances are wanted at explosion, end of
network integration, grain-formation time, or after decay. A final-file name
does not answer that question.

Before processing every particle, resolve model association and trajectory
physics, choose a common epoch or explicitly retain differing epochs, establish
archive integrity and completed-job coverage, then stream the data into
partitioned outputs in user scratch. The current bounded audit is a discovery
and validation product; it does not claim every archive byte or every PRISM
particle has been audited.
