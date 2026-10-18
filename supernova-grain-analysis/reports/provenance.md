# Hydro provenance and adversarial checks

Read-only inspection on 2026-09-25. These findings concern simulation/SPH
particles, not condensed dust grains. Paths below are source data and must
remain read-only. This audit inspected headers, source code and bounded samples;
it does not establish complete equivalence of copied datasets.

## Raw hydro locations

| Requested model | Observed source directory | Numeric-suffix binary snapshots |
|---|---|---|
| 50AM (directory spelling `50Am`) | `/data/pyoung3/supernova_data/50Am/run3g_50Am6` | 61: `run3g_50Am6_sph.0000` through `.2000` |
| jet3b | `/data/pyoung3/supernova_data/jet3b/jet3b` | 61: `run3g_50Am_jet3b_sph.0000` through `.8000` |
| cco2 | `/data/pyoung3/gsvance/supernova_data/cco2/cco2` | 190: `r3g_1M_cco2_sph.1000` through `.19900` |
| vconvL | `/data/pyoung3/gsvance/supernova_data/vconvL/vconvL_snsph` | No numeric-suffix binary snapshots in this directory; text `.out` files exist |
| Related vconv hydro | `/data/pyoung3/gsvance/supernova_data/vconv/vconv_snsph` | 170: `run3g_50_1milvc2_sph.0000` through `.3380` |

A later binary snapshot is outside the last directory:
`/data/pyoung3/gsvance/supernova_data/vconv/run3g_50_1milvc2_sph.4200`.
Its header has `npart=945342`, `tpos=12283.2`, `timeCF=100`.

The first particle, ID 0, of `vconvL_snsph/...0000.out` agrees, at printed
precision, with the related vconv `.0000` binary in position, temperature,
density, velocity and mass. The first particle of `vconvL_snsph/...4200.out`
likewise agrees with the separate vconv `.4200` binary. This supports shared
hydro provenance for those samples. It does **not** establish that `vconv`,
`vconvL`, and collaborator PRISM results are interchangeable. The vconvL
`link_burn_query_output.bash` links legacy vconv burn-query logs; their analysis
tables are separate and have different sizes. Preserve the original model label
and path until a complete particle/time comparison establishes an association.

## File meanings and fields

An SDF snapshot is a text header followed by binary particle records. Its suffix
is an iteration identifier. Its header supplies `tpos`, particle count, field
types, and sometimes conversion factors. Do not assign a universal binary dtype:
50Am/jet3b headers have particle-level `p1..p22`, `m1..m22`, while sampled
cco2/vconv files declare header-level `p1..p20`, `n1..n20` and have 200-byte
particle records. Derive isotope identity from the appropriate declaration.

The observed 14-column, headerless VRout mapping is:

| Column (one based) | Meaning |
|---|---|
| 1 | Particle `ident` |
| 2–4 | `x,y,z` |
| 5 | `temp` |
| 6 | `rho` |
| 7 | `pr` (pressure) |
| 8–10 | `vx,vy,vz` |
| 11 | Particle mass |
| 12 | `Y_el` |
| 13 | Position magnitude `sqrt(x*x+y*y+z*z)` |
| 14 | **Speed** `sqrt(vx*vx+vy*vy+vz*vz)`, despite the variable name `vr` |

Evidence: `/data/pyoung3/gsvance/data_mining/cco2-SDF-reader.c:124–129`
and `SDF-reader.c:104–109`. Actual first records of cco2 `.1000` and vconv
`.0000` match their VRout text at printed precision. No unit conversion is
applied in these exporters. Column 14 is not signed radial velocity; calculate
the latter as `(x*vx+y*vy+z*vz)/r`, with an explicit policy for `r=0`.

The vconvL `.out` sample has a different, explicit 17-column header:
`ID, X_Pos, Y_Pos, Z_Pos, Temp, U, U_dot, rho, V_x, V_y, V_z, A_x, A_y, A_z, h, Mass, Y_e`.
Here density is column **8**, not column 6. Filename extension alone cannot
choose the parser. Existing `.out` files may represent selections rather than
every snapshot/particle; compare IDs, not row positions.

## Unit evidence and concrete discrepancies

`/data/pyoung3/gsvance/data_mining/sn_utils.py:56–63` records SNSPH conversion
factors from an initial control file:

- mass: `1.9889e27 g` per code mass;
- length: `6.955e10 cm` per code length;
- time: `100 s` per code time;
- density: `massCF / lengthCF**3 = 5.911824660348532e-6 g/cm^3` per code density;
- velocity: `lengthCF / timeCF = 6.955e8 cm/s` per code velocity.

Sampled cco2 headers explicitly contain `massCF`, `lengthCF`, `timeCF`.
Sampled vconv headers instead spell the length field **`lenCF`**. Readers must
recognize both, reject inconsistent duplicate declarations, and record whether
any factor was inferred. Sampled 50Am/jet3b headers omit these factors; the
shared code provides evidence for defaults, not a header-level guarantee.
`/data/pyoung3/snsphtools/SDFtoASCII.c:156–158` reads header factors or falls
back to these same values. `maketraj.c:163–164,243–247` uses `tpos*100`, the
same length conversion, and the same density conversion.

The collaborator scripts have material discrepancies that must be resolved
before using their trajectory products quantitatively:

1. `/scratch/svierlin/prism-v02/generate_trajectories_jet3b.sh:63–66` writes
   the snapshot suffix as seconds. Raw jet3b `.8000` instead has
   `iter=8000`, `tpos=964535`: using the documented default gives
   `96,453,500 s`, not `8000 s`.
2. `generate_trajectories_cco2.sh:68` onward uses a hardcoded time array
   beginning `1.47e-2`, and lines 279–280 format it directly as time. The raw
   `.1000` header gives `tpos=.014655`, `timeCF=100`: physical time is
   `1.4655 s`. The list begins at rounded **code time**, not physical seconds.
3. `generate_trajectories_vconv.sh:68` onward similarly begins `0,3.5e-6,...`.
   The `.0020` header explicitly has `tpos=3.5e-6`, `timeCF=100`, hence
   `0.00035 s`. Lines 259–260 use the unconverted value and round to three
   significant figures.
4. The jet3b script line 97, cco2 line 312 and vconv line 292 multiply raw
   density by `1.42e-6`. The header-supported factor is **4.163256803 times
   larger**. The origin of `1.42e-6` remains unexplained. Do not silently
   reproduce it or describe these products as validated CGS trajectories.
5. `extract_ids_vconvL.sh:25,48–59` scans vconvL `.out` in increments of 100,
   although the directory contains early files in increments of 20. It also
   converts every first-column value numerically without skipping the header,
   which turns `ID` into 0. Its shared `input/conditions/id_list.txt` filename
   does not identify a model. These facts do not prove a wrong ID list was
   actually used, but its provenance must be checked.

## Remaining checks

- Validate vconvL association beyond sampled endpoints, including selected
  particle populations and downstream abundance products.
- Verify the original control files for the 50Am/jet3b default unit factors.
- Determine which collaborator products were generated by which script
  revision; script inspection alone does not date or validate an output.
- Preserve explicitly whether abundances are initial, hydro-network,
  postprocessed, or decay-adjusted. Related filenames do not establish that
  distinction.
