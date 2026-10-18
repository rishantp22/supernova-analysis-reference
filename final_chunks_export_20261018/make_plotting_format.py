"""Create the requested comma-space table, preserving original data files."""
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
SOURCE = ROOT / 'final_chunks.out'
REFERENCE = Path('/data/pyoung3/supernova_data/50Am/analysis/50Am_plotting.out')
OUTPUT = ROOT / '50Am_plotting_all_isotopes.out'
NAMES = (ROOT / 'List of isotopes_elements_variables.txt').read_text().splitlines()
PHYSICAL = NAMES[:16]
REFERENCE_NAMES = {'Id': 'id', 'Density': 'density', 'Peak temp': 'peak temp',
                   'Peak density': 'peak density'}
FACTORS = {**dict.fromkeys(['x','y','z','h'], 6.955e10),
           **dict.fromkeys(['vx','vy','vz'], 6.955e8),
           **dict.fromkeys(['ax','ay','az'], 6.955e6),
           'mass': 1.9889e27, 'Density': 1.9889e27 / (6.955e10 ** 3)}
UNITS = {**dict.fromkeys(['x','y','z','h'], 'cm'),
         **dict.fromkeys(['vx','vy','vz'], 'cm/s'),
         **dict.fromkeys(['ax','ay','az'], 'cm/s^2'),
         'mass':'g', 'Density':'g/cm^3', 'Peak temp':'K', 'Peak density':'g/cm^3',
         'Y_{e}':'dimensionless; progenitor electron fraction as in reference table', 'Id':'ID'}


def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True)


def stat(path):
    s = path.stat()
    return {'path':str(path), 'bytes':s.st_size, 'mtime_ns':s.st_mtime_ns}


def output_name(name):
    if name in PHYSICAL:
        return REFERENCE_NAMES.get(name, name)
    return name if name == 'X_{n}' else 'X_{' + name + '}'


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    WORK.mkdir(exist_ok=False)
    before = [stat(p) for p in [SOURCE, REFERENCE, ROOT / 'List of isotopes_elements_variables.txt']]
    original_validation = json.loads((WORK.parent / 'validation.json').read_text())
    log('Reading the previously validated 298-column table')
    data = pd.read_csv(SOURCE, sep='\t', dtype={'Id':'int64'})
    assert list(data) == NAMES
    assert len(data) == original_validation['rows'] and data.Id.is_unique
    data.sort_values('Id', inplace=True)
    data.set_index('Id', drop=False, inplace=True)
    ref_names = [REFERENCE_NAMES.get(n,n) for n in PHYSICAL]
    probe_names = ['X_{4He}', 'X_{12C}', 'X_{16O}', 'X_{56Fe}']
    log('Reading reference physical fields and checking the abundance-source distinction')
    reference = pd.read_csv(REFERENCE, skipinitialspace=True, usecols=ref_names + probe_names,
                            dtype={**dict.fromkeys(ref_names + probe_names, str), 'id':'int64'},
                            keep_default_na=False).set_index('id', drop=False)
    assert reference.index.is_unique
    assert data.index.isin(reference.index).all()
    reference_rows = len(reference)
    reference = reference.reindex(data.index)
    physical_checks = {}
    for name, factor in FACTORS.items():
        values = pd.to_numeric(reference[REFERENCE_NAMES.get(name,name)]).to_numpy()
        expected = data[name].to_numpy() * factor
        assert np.isfinite(values).all()
        # Reference printing can introduce slight decimal-rounding differences.
        ok = np.isclose(values, expected, rtol=1e-9, atol=1e-300)
        physical_checks[name] = int((~ok).sum())
        if not ok.all(): raise ValueError(f'Reference physical mismatch: {name}, {physical_checks[name]} rows')
    abundance_checks = {}
    for reference_name in probe_names:
        isotope = reference_name[3:-1]
        old = pd.to_numeric(reference[reference_name]).to_numpy()
        new = data[isotope].to_numpy()
        match = np.isclose(old, new, rtol=1e-5, atol=1e-30)
        abundance_checks[isotope] = {'matching_rows':int(match.sum()), 'different_rows':int((~match).sum())}
    reference.drop(columns=probe_names, inplace=True)
    # Physical definitions deliberately follow the user's actual reference file.
    # Its Y_e is progenitor Y_e, and its peaks use the legacy early-snapshot set.
    semantic_changes = {}
    for name, factor in [('Peak temp',1.),('Peak density',FACTORS['Density']),('Y_{e}',1.)]:
        actual = pd.to_numeric(reference[REFERENCE_NAMES.get(name,name)]).to_numpy()
        expected = data[name].to_numpy() * factor
        semantic_changes[name] = int((~np.isclose(actual, expected, rtol=1e-5, atol=1e-300)).sum())
    header = [output_name(n) for n in NAMES]
    assert len(header) == 298 and len(set(header)) == 298
    missing = {output_name(n):int(data[n].isna().sum()) for n in NAMES[16:]}
    column_map = []
    for name in NAMES:
        if name in PHYSICAL:
            source = str(REFERENCE)
            definition = 'Copied verbatim by matching id from the reference plotting table'
        else:
            source = str(SOURCE)
            definition = ('Unavailable in the supplied chunks; NA retained'
                          if data[name].isna().all() else 'Unchanged last-category chunk isotope mass fraction')
        column_map.append({'requested_name':name, 'output_name':output_name(name),
                           'units':UNITS.get(name,'dimensionless mass fraction'),
                           'source':source, 'definition':definition})
    pd.DataFrame(column_map).to_csv(WORK / 'column_map.csv', index=False, mode='x')
    partial = WORK / 'table.partial'
    digest = hashlib.sha256()
    with partial.open('xb') as stream:
        first_line = (', '.join(header) + '\n').encode()
        stream.write(first_line); digest.update(first_line)
        for start in range(0, len(data), 10000):
            block = data.iloc[start:start+10000].copy()
            for name in PHYSICAL:
                if name != 'Id': block[name] = reference[REFERENCE_NAMES.get(name,name)].iloc[start:start+len(block)].to_numpy()
            block.columns = header
            payload = block.to_csv(index=False, header=False, na_rep='NA', float_format='%.15g').replace(',', ', ').encode()
            stream.write(payload); digest.update(payload)
            if start % 100000 == 0: log(f'Writing rows: {start+len(block):,}/{len(data):,}')
    log('Verifying every written row, header, ID and column count')
    count = 0
    expected_ids = data.index.to_numpy()
    with partial.open() as stream:
        actual_header = stream.readline().rstrip('\n').split(', ')
        assert actual_header == header
        for line in stream:
            values = line.rstrip('\n').split(', ')
            assert len(values) == len(header)
            ident = int(values[0])
            assert ident == expected_ids[count]
            count += 1
    assert count == len(data)
    after = [stat(p) for p in [SOURCE, REFERENCE, ROOT / 'List of isotopes_elements_variables.txt']]
    assert before == after
    os.link(partial, OUTPUT)
    partial.unlink()
    report = {'output':str(OUTPUT), 'rows':count, 'columns':len(header),
              'delimiter':'comma followed by a space', 'ID_order':'ascending numeric',
              'output_bytes':OUTPUT.stat().st_size, 'output_sha256':digest.hexdigest(),
              'requested_columns_mapped_in_exact_original_order':True,
              'physical_units':'CGS', 'reference_total_IDs':reference_rows,
              'reference_IDs_absent_from_chunks':reference_rows-count,
              'physical_conversion_mismatches':physical_checks,
              'composition_category':'last', 'composition_source':'validated supplied AM chunks',
              'reference_vs_chunks_abundance_comparison':abundance_checks,
              'reference_vs_previous_export_semantic_differences':semantic_changes,
              'missing_cells_by_composition_column':missing,
              'input_metadata_unchanged':True, 'sources':before,
              'prior_export_sha256':original_validation['output_sha256']}
    with (WORK / 'validation.json').open('x') as stream:
        json.dump(report, stream, indent=2)
    notes = f'''50Am_plotting_all_isotopes.out

{count:,} unique particle IDs; 298 columns; one header line; comma-space separated.
Rows are sorted by numeric ID. The original requested list determines column order.
Names follow the example: id, density, peak temp, peak density, X_{{H}}, X_{{1H}}, etc.
column_map.csv maps every original requested name to the output name and source.

The first 16 physical fields are copied by ID from the actual reference file:
{REFERENCE}
Positions and h: cm; velocities: cm/s; acceleration: cm/s^2; mass: g;
density and peak density: g/cm^3; peak temp: K. Y_{{e}} is the PROGENITOR
electron fraction used by the reference generator, not final snapshot Y_e.
Peaks follow the reference's early-snapshot definition and printed precision.
The length/mass/density conversions agree with the validated scratch physical
file for all matched IDs. Local SNSPH conversion conventions are preserved.

All 222 requested isotope columns supplied by the chunks retain their last-category
C_sum values. No isotope value was replaced by a legacy burn-query abundance.
The reference plotting file is an older burn-query product; its composition is
different from the supplied PRISM chunks. Mixing its elemental totals with the
new isotope values would not form one consistent composition.
Consequently all 59 elemental-total columns and X_{{n}} remain NA, as their
complete values cannot be recovered from the selected-isotope chunks. NA is
missing information, never zero. Existing source zeros are preserved.

The reference has {reference_rows:,} IDs, of which {reference_rows-count:,} have no chunk abundances.
The output contains every one of the {count:,} chunk IDs. Example IDs 0, 2 and 3
are among those missing from the chunks, so they cannot be extended with the
requested chunk isotopes. The prior hydro_ids_without_chunks.csv lists those IDs.

The isotope sum anomaly at ID 7753 remains unchanged. No renormalization or decay
correction is applied. Abundance epochs are retained in the prior grain_epochs.csv.
All source files, the original plotting table and final_chunks.out remain unchanged.
'''
    with (WORK / 'README.txt').open('x') as stream: stream.write(notes)
    with (WORK / 'preview.out').open('x') as preview, OUTPUT.open() as full:
        for _ in range(6): preview.write(full.readline())
    log(f'COMPLETE: {OUTPUT}, {count:,} rows, 298 columns')
    print(json.dumps({k:v for k,v in report.items() if k != 'missing_cells_by_composition_column'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
