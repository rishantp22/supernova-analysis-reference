# Complete chunk integrity audit

All 938 A/B chunk pairs were inspected, not sampled. A files total 31,749,696,840 bytes and 624,194,514 records; B files total 19,098,031,272 bytes and the same record count. The complete grid is 937,229 unique IDs × three categories × 222 species. Categories are nearest saved time to 100 s, nearest saved time to 1000 s, and last saved time; they are not guaranteed exact target epochs.

There are zero duplicate cross-chunk IDs, duplicate keys, missing ID/category/species cells, inconsistent A/B keys, negative/nonfinite abundances, or B mass-weighting discrepancies at relative tolerance 2e-5. A's column named `A` is Z and column `B` is nuclear mass number. B abundances use native particle-mass units of 1e-6 solar masses, not grams.

Every valid new plotting value exactly equals its last-category A abundance. The two chunk IDs absent from valid new-file IDs are 23 and 193; the latter has malformed token `y193` in the plotting file. No automatic repair was applied. 8,399 valid old IDs are absent from the chunks; see the stage ledger for differences against candidate inputs and independent EOC coverage.

Last saved times range from 0.936977 to 166000 seconds, with median 166000. The maximum last-category selected sum is 1.977369060913544 for ID 7753. There is no model/run identifier in a chunk row; complete grids alone cannot certify run identity or successful nuclear integration.

Independent archived final EOC comparisons found 8,456 finite mismatches in 176 particles using rtol=5.1e-6 and atol=0. All these particles have chunk endpoints at or before 2000 s. Particle 7753 has approximately doubled values for 48 finite discrepant species. A/B internal consistency is therefore compatible with an upstream extraction error. The producer's printed-time summation is a demonstrated mechanism capable of merging repeated blocks, but the actual affected 50AM time series were not recovered, so this mechanism remains an evidence-supported hypothesis.

During the first pass, 24 chunks changed filesystem metadata while being read. Their failed attempts are preserved. The retry required matching content hashes when metadata changed; all 938 pairs completed successfully. End-of-audit source stability checks rehashed changed-metadata sources and rehashed both plotting files. This is not evidence that the initial metadata changes altered abundance bytes. Receipts include paths, sizes, metadata and SHA-256.

Detailed machine-readable evidence: `tables/chunk_integrity.csv`, `tables/all_chunk_particles.npz`, `tables/eoc_mismatches_with_chunk_epoch.csv`, `tables/eoc_plot_mismatch_cells.csv`, `provenance/chunks/`, `provenance/source_stability.json`.
