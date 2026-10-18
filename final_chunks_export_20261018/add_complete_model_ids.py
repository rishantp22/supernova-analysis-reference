"""Include every reference-model ID without inventing missing compositions."""
import csv
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd

ROOT = Path('/scratch/rpraka15')
WORK = Path(__file__).resolve().parent / 'plotting_format'
MATCHED = ROOT / '50Am_plotting_all_isotopes.out'
REFERENCE = Path('/data/pyoung3/supernova_data/50Am/analysis/50Am_plotting.out')
OUTPUT = ROOT / '50Am_plotting_complete.out'


def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True)


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    prior = json.loads((WORK / 'validation.json').read_text())
    with MATCHED.open() as stream: header = stream.readline().rstrip('\n').split(', ')
    assert len(header) == 298
    missing = set(pd.read_csv(WORK.parent / 'hydro_ids_without_chunks.csv').ID)
    missing_rows, all_ids = {}, []
    ref_stat = REFERENCE.stat()
    log(f'Reading physical fields for {len(missing):,} additional model IDs')
    for block in pd.read_csv(REFERENCE, skipinitialspace=True, usecols=header[:16],
                             dtype=str, keep_default_na=False, chunksize=50000):
        ids = pd.to_numeric(block.id).to_numpy(dtype=np.int64)
        all_ids.append(ids)
        selected = block.loc[np.isin(ids, list(missing)), header[:16]]
        for row in selected.itertuples(index=False, name=None):
            ident = int(row[0])
            assert ident not in missing_rows
            missing_rows[ident] = (', '.join(list(row) + ['NA'] * (len(header)-16)) + '\n').encode()
    all_ids = np.concatenate(all_ids)
    assert len(np.unique(all_ids)) == len(all_ids)
    assert set(missing_rows) == missing
    ordered_ids = np.sort(all_ids)
    missing_order = sorted(missing)
    partial = WORK / 'complete_model.partial'
    digest, input_digest = hashlib.sha256(), hashlib.sha256()
    count, at = 0, 0
    coverage_path = WORK / 'complete_model_coverage.csv'
    with MATCHED.open('rb') as source, partial.open('xb') as out, coverage_path.open('x') as coverage:
        coverage.write('id,chunk_composition_available\n')
        first = source.readline(); input_digest.update(first)
        out.write(first); digest.update(first)
        for line in source:
            input_digest.update(line)
            ident = int(line.partition(b', ')[0])
            while at < len(missing_order) and missing_order[at] < ident:
                other = missing_order[at]
                assert other == ordered_ids[count]
                payload = missing_rows[other]
                out.write(payload); digest.update(payload)
                coverage.write(f'{other},false\n')
                count += 1; at += 1
            assert ident == ordered_ids[count] and ident not in missing
            out.write(line); digest.update(line)
            coverage.write(f'{ident},true\n')
            count += 1
        while at < len(missing_order):
            ident = missing_order[at]
            assert ident == ordered_ids[count]
            payload = missing_rows[ident]
            out.write(payload); digest.update(payload)
            coverage.write(f'{ident},false\n')
            count += 1; at += 1
    assert input_digest.hexdigest() == prior['output_sha256']
    assert count == len(all_ids) == prior['reference_total_IDs']
    assert count-len(missing) == prior['rows']
    assert (REFERENCE.stat().st_size, REFERENCE.stat().st_mtime_ns) == (ref_stat.st_size, ref_stat.st_mtime_ns)
    log('Verifying all complete-model rows and the example IDs')
    check_digest = hashlib.sha256()
    seen_count = 0
    examples = {}
    with partial.open('rb') as stream:
        first = stream.readline(); check_digest.update(first)
        assert first.decode().rstrip('\n').split(', ') == header
        for line in stream:
            check_digest.update(line)
            assert line.count(b', ') == 297
            ident = int(line.partition(b', ')[0])
            assert ident == ordered_ids[seen_count]
            if ident in [0, 2, 3]:
                values = line.decode().rstrip('\n').split(', ')
                assert all(x == 'NA' for x in values[16:])
                examples[ident] = dict(zip(header[:16], values[:16]))
            seen_count += 1
    assert seen_count == count and check_digest.hexdigest() == digest.hexdigest()
    os.link(partial, OUTPUT)
    partial.unlink()
    report = {'output':str(OUTPUT), 'rows':count, 'columns':len(header),
              'all_reference_model_IDs_included':True, 'duplicate_IDs':0,
              'rows_with_chunk_composition':prior['rows'], 'rows_without_chunk_composition':len(missing),
              'unavailable_composition_columns_for_all_rows':60,
              'unavailable_composition_columns_for_missing_chunk_rows':282,
              'delimiter':'comma-space', 'physical_units':'CGS', 'sort_order':'numeric id ascending',
              'original_matched_rows_preserved_byte_for_byte':True,
              'output_sha256':digest.hexdigest(), 'output_bytes':OUTPUT.stat().st_size,
              'example_physical_rows_verified':examples,
              'source_files_unchanged':True}
    with (WORK / 'complete_model_validation.json').open('x') as stream:
        json.dump(report, stream, indent=2)
    with (WORK / 'complete_model_preview.out').open('x') as preview, OUTPUT.open() as source:
        for _ in range(6): preview.write(source.readline())
    notes = f'''50Am_plotting_complete.out

All {count:,} AM IDs from the actual 50Am_plotting.out reference are included,
sorted numerically, with 298 comma-space-separated columns. Header spelling follows
the example (id, density, X_{{4He}}, etc.); column order follows the requested list.

The first 16 fields are copied by ID from the reference plotting table, in CGS units.
Y_{{e}} is its progenitor electron fraction. Peaks retain its early-snapshot definition.
{prior['rows']:,} IDs have all 222 requested isotope values supplied by the last-category
chunks. Those values are unchanged from the validated chunk export.
{len(missing):,} IDs, including 0, 2 and 3, have no supplied chunk compositions: their
physical fields are present and all 282 composition fields are explicitly NA.
complete_model_coverage.csv identifies these rows.

For all rows, 59 full elemental totals and X_{{n}} remain NA because the selected
isotope chunks cannot determine them. NA is unavailable information, not zero.
The older reference's composition belongs to a different abundance calculation;
it was not mixed into the new chunk isotope values. No missing data were invented.
All 298 requested variables are represented, but a fully populated numerical table
cannot be obtained from these chunks alone.

Original sources and earlier outputs remain unchanged. See column_map.csv for the
field mapping and README.txt for physical definitions and source distinctions.
complete_model_validation.json records counts, checks and SHA-256.
'''
    with (WORK / 'complete_model_README.txt').open('x') as stream: stream.write(notes)
    log(f'COMPLETE: {OUTPUT}; {count:,} rows, 298 columns')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
