# Particle-output coverage and archive limitations

The output census is a list of direct entry names, not proof that every particle
finished post-processing. Numeric names are candidate SPH particle IDs. The
census deliberately avoids recursively listing hundreds of thousands of
particle directories or issuing one metadata request per particle.

## Source layout

| Requested model | Observed output root | Representation |
| --- | --- | --- |
| 50AM | `/data/pyoung3/groupdata/svierlin/50Am_data` | Ten `eoc_batch_*.tar.gz` archives |
| jet3b | `/data/pyoung3/groupdata/svierlin/jet3b_data/jet3b_outputs` | Numeric direct entry names |
| cco2 | `/data/pyoung3/groupdata/svierlin/cco2_data/cco2_outputs` | Numeric direct entry names |
| vconvl (unconfirmed association with observed `vconv`) | `/data/pyoung3/groupdata/svierlin/vconv_data/vconv_outputs` | Numeric direct entry names |

`vconvl` is the requested label; `vconv` is the observed folder spelling. Their
equivalence requires evidence from configuration or provenance, beyond the
folder-name association. Likewise, the observed archive folder spells `50Am`.

## Reproduce the census

Run `python scripts/inventory.py` from the repository. It writes only to
`/scratch/rpraka15/supernova-grain-analysis/inventory` by default:

- `inventory.json`: counts, bounds, nonnumeric names, access errors, archive
  sizes, and bounded log findings. A failed census is explicitly marked partial.
- `jet3b_particle_ids.csv`, `cco2_particle_ids.csv`, and
  `vconv_particle_ids.csv`: exact observed numeric names and source paths.
- `sample_ids.json`: up to 24 deterministic IDs per model, selected at evenly
  spaced positions in the numerically sorted ID list. This is a structural
  debugging sample, not a random or mass-weighted scientific sample.
- `archive_log_excerpts.txt`: first and last 4096 bytes of at most the first 96
  sorted entries in each of four designated log folders.

Directory enumeration is not an atomic snapshot. Files may change during the
audit, and counts describe the observation time recorded in the inventory.
Numeric entries have not all been checked to be directories, to contain valid
output, or to represent unique physical particles across different models.

An early intermediate census file named `vconvl_particle_ids.csv` is preserved
under the owner's no-deletion instruction. Use the corrected
`vconv_particle_ids.csv` and current inventory JSON, which retain the observed
label. That old filename does not establish a vconvL association.

The names-only census on 2026-09-25 observed 4,250 numeric entries for `jet3b`,
137,061 for `cco2`, and 167,217 for observed `vconv`, with no census access
errors. These are output-root entry counts, not completeness fractions: the
denominators of intended runs and the definition of successful output still
need independent verification. Machine-readable records use the observed
`vconv` label rather than silently equating it with requested `vconvl`.

## Archive health

The initial listing observed ten archives totaling 2,167,464,837,911 bytes
(approximately 2.17 TB in decimal units), including three with sizes exactly
divisible by 512 bytes: batch 0001 (524,214,599,680 bytes), batch 0008
(376,129,716,224 bytes), and batch 0009 (385,585,250,304 bytes). This is a
diagnostic clue only: file size divisibility does not establish truncation or
validity. The archives are too large to casually decompress in full on a login
node. Archive records therefore retain `integrity_verified: false`.

Bounded log reads found cancellations on May 12 for jobs 52679133, 52679134,
52679140, 52679141, and 52679142. Jobs 52679135 through 52679139 logged finish
messages. On May 13, rerun jobs 52754098 and 52754100 logged finish messages,
while jobs 52754099, 52754101, and 52754102 reported `gzip: stdout: No space
left on device` and tar errors. Extraction job 62563401 also logged a failure.
These are positive evidence of interrupted archive operations, not merely
suspicious file sizes. Paired `.out` logs identify the rerun targets:

| Job | Batch | Intended subfolders in log | Observed outcome |
| --- | --- | ---: | --- |
| 52754098 | 0010 | 37,229 | Finish message |
| 52754099 | 0001 | 100,000 | No-space gzip/tar error; no finish message |
| 52754100 | 0007 | 100,000 | Finish message |
| 52754101 | 0008 | 100,000 | No-space gzip/tar error; no finish message |
| 52754102 | 0009 | 100,000 | No-space gzip/tar error; no finish message |

The failing batches match the three archive sizes divisible by 512, which
strengthens the concern that batches 0001, 0008, and 0009 are incomplete. The
current files have not received full integrity checks; later replacement or
partial recoverability cannot be resolved from these sampled logs.

The extraction error log reports available space of 939,772,608,512 bytes
against an estimated requirement of 1,572,643,799,040 bytes, followed by pigz
reporting empty input and tar reporting that its input was not a tar archive.
This proves that extraction attempt failed; it does not by itself identify an
additional corrupt archive. Exact log paths and sampled lines are retained in
the scratch inventory and excerpt file.

Successfully reading an archive prefix can identify an early member and its
format; it cannot verify the gzip trailer, tar terminator, remaining members,
or whether all intended particles were archived. Logs may provide independent
evidence of interrupted jobs, but do not by themselves establish which archive
members are usable. The original manifest, intended batch membership, job
completion state, and full integrity verification remain separate questions.

The log sampler does not read beyond the first/last 4096 bytes of each chosen
file. Error-looking lines are search hits requiring interpretation, and absence
of such a line in a sample is not proof of successful processing. All source
directories are read-only for this project; nothing is extracted there.
