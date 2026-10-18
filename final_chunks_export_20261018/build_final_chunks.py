#!/usr/bin/env python3
"""Export existing AM chunks by explicit ID/species joins; never edit sources."""
import argparse
import concurrent.futures as futures
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path('/scratch/rpraka15')
WORK = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'supernova-grain-analysis/code'))
from readers import open_sdf, request_kind

CATEGORIES = ['closest_to_100', 'closest_to_1000', 'last']
CAT_TYPE = pd.CategoricalDtype(CATEGORIES)
PHYSICAL = {'Id': 'ID', 'x': 'X_Pos', 'y': 'Y_Pos', 'z': 'Z_Pos',
            'vx': 'V_x', 'vy': 'V_y', 'vz': 'V_z', 'ax': 'A_x',
            'ay': 'A_y', 'az': 'A_z', 'mass': 'Mass', 'h': 'h',
            'Density': 'rho', 'Y_{e}': 'Y_e'}
RAW_FIELDS = {'x': 'x', 'y': 'y', 'z': 'z', 'vx': 'vx', 'vy': 'vy',
              'vz': 'vz', 'ax': 'ax', 'ay': 'ay', 'az': 'az',
              'mass': 'mass', 'h': 'h', 'Density': 'rho', 'Y_{e}': 'Y_el'}


def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True)


def fingerprint(path):
    stat = path.stat()
    return {'path': str(path), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def load_csv(path, **kwargs):
    before = fingerprint(path)
    data = path.read_bytes()
    frame = pd.read_csv(io.BytesIO(data), **kwargs)
    if fingerprint(path) != before:
        raise ValueError(f'Source changed during read: {path}')
    before['sha256'] = hashlib.sha256(data).hexdigest()
    return frame, before


def save_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def get_physical():
    """Keep text-export values exactly and verify their AM association by ID."""
    source = ROOT / 'run3g_50Am6_sph.2000.out'
    frame, source_info = load_csv(source, dtype=str, skipinitialspace=True,
                                  keep_default_na=False)
    if set(PHYSICAL.values()) - set(frame.columns):
        raise ValueError('Physical source lacks required columns')
    ids = pd.to_numeric(frame.ID, errors='raise').to_numpy(dtype=np.int64)
    if len(np.unique(ids)) != len(ids):
        raise ValueError('Duplicate physical IDs')
    frame.index = ids
    raw_path = Path('/data/pyoung3/supernova_data/50Am/run3g_50Am6/run3g_50Am6_sph.2000')
    info, raw = open_sdf(raw_path)
    raw_ids = np.asarray(raw['ident'], dtype=np.int64)
    if len(np.unique(raw_ids)) != len(raw_ids):
        raise ValueError('Duplicate raw snapshot IDs')
    order = np.argsort(raw_ids)
    ordered = raw_ids[order]
    loc = np.searchsorted(ordered, ids)
    if (loc >= len(ordered)).any() or not np.array_equal(ordered[loc], ids):
        raise ValueError('Physical text contains IDs absent from raw AM snapshot')
    raw_index = order[loc]
    check = {}
    for output, raw_name in RAW_FIELDS.items():
        text_values = pd.to_numeric(frame[PHYSICAL[output]], errors='raise').to_numpy()
        raw_values = np.asarray(raw[raw_name][raw_index], dtype=np.float64)
        good = np.isfinite(text_values) & np.isfinite(raw_values)
        good &= np.isclose(text_values, raw_values, rtol=5.1e-6, atol=1e-300)
        check[output] = int((~good).sum())
        if not good.all():
            raise ValueError(f'Physical text/raw mismatch in {output}: {check[output]} rows')
    verification = {'rows': len(frame), 'raw_snapshot': fingerprint(raw_path),
                    'raw_rows': len(raw), 'raw_iteration': int(info['meta']['iter']),
                    'raw_time_code': info['meta']['tpos'],
                    'raw_text_mismatches': check,
                    'comparison_rtol': 5.1e-6, 'source': source_info}
    log(f'Physical data: {len(frame):,} unique IDs; all 13 numeric fields agree with raw AM at text precision')
    return frame, verification


def get_peaks():
    """Maxima over every available saved AM snapshot, using explicit IDs."""
    folder = Path('/data/pyoung3/supernova_data/50Am/run3g_50Am6')
    paths = sorted(p for p in folder.glob('run3g_50Am6_sph.*')
                   if re.fullmatch(r'run3g_50Am6_sph\.\d+', p.name))
    if not paths:
        raise ValueError('No AM snapshots available for peak fields')
    cache = WORK / 'peaks.npz'
    metadata_path = WORK / 'peak_sources.json'
    if cache.exists() and metadata_path.exists():
        meta = json.loads(metadata_path.read_text())
        if [fingerprint(p) for p in paths] != [r['source'] for r in meta['snapshots']]:
            raise ValueError('Snapshot sources changed since peak cache was built')
        with np.load(cache) as data:
            result = pd.DataFrame({name: data[name] for name in data.files if name != 'Id'},
                                  index=data['Id'])
        return result, meta
    max_t = np.full(1000000, -np.inf)
    max_rho = np.full(1000000, -np.inf)
    counts = np.zeros(1000000, dtype=np.int32)
    records = []
    for n, path in enumerate(paths, 1):
        before = fingerprint(path)
        info, raw = open_sdf(path)
        ids = np.asarray(raw['ident'], dtype=np.int64)
        if len(np.unique(ids)) != len(ids) or ids.min() < 0:
            raise ValueError(f'Duplicate/invalid IDs in {path}')
        if ids.max() >= len(max_t):
            extra = int(ids.max()) + 1 - len(max_t)
            max_t = np.pad(max_t, (0, extra), constant_values=-np.inf)
            max_rho = np.pad(max_rho, (0, extra), constant_values=-np.inf)
            counts = np.pad(counts, (0, extra))
        t = np.asarray(raw['temp'], dtype=np.float64)
        rho = np.asarray(raw['rho'], dtype=np.float64)
        if not np.isfinite(t).all() or not np.isfinite(rho).all() or (t < 0).any() or (rho < 0).any():
            raise ValueError(f'Invalid temperature/density in {path}')
        max_t[ids] = np.maximum(max_t[ids], t)
        max_rho[ids] = np.maximum(max_rho[ids], rho)
        counts[ids] += 1
        if fingerprint(path) != before:
            raise ValueError(f'Snapshot changed during read: {path}')
        records.append({'source': before, 'iteration': int(info['meta']['iter']),
                        'time_code': info['meta']['tpos'], 'rows': len(ids),
                        'header_sha256': hashlib.sha256(info['header'].encode()).hexdigest()})
        del raw
        log(f'Peak scan: {n}/{len(paths)} snapshots')
    ids = np.flatnonzero(counts)
    result = pd.DataFrame({'Peak temp': max_t[ids], 'Peak density': max_rho[ids],
                           'hydro_snapshot_count': counts[ids]}, index=ids)
    meta = {'definition': 'Independent maxima over all available saved AM hydro snapshots; not a continuous-time peak or PRISM input-trajectory peak',
            'temperature_units': 'K, as stored in SDF temp',
            'density_units': 'native code density, same as output Density',
            'snapshots': records}
    with cache.open('xb') as handle:
        np.savez(handle, Id=ids, **{name: result[name].to_numpy() for name in result})
    save_json(metadata_path, meta)
    return result, meta


def encoded_keys(frame):
    return ((frame.ID.to_numpy(dtype=np.int64) * 3
             + frame.category.cat.codes.to_numpy(dtype=np.int64)) * 65536
            + frame.A.to_numpy(dtype=np.int64) * 512 + frame.B.to_numpy(dtype=np.int64))


def process_chunk(n, args, columns, physical, peaks):
    path_a = ROOT / 'final_chunks' / f'chunk_{n}.csv'
    path_b = ROOT / 'final_chunks/final_chunks_B' / f'chunk_{n}_B.csv'
    common = {'ID': 'int64', 'category': CAT_TYPE, 'A': 'int16', 'B': 'int16'}
    a, meta_a = load_csv(path_a, dtype={**common, 'actual_timestep': 'int64',
                                      'actual_time': 'float64', 'C_sum': 'float64'})
    b, meta_b = load_csv(path_b, dtype={**common, 'mass_times_Csum': 'float64'})
    if list(a.columns) != ['ID','category','actual_timestep','actual_time','A','B','C_sum']:
        raise ValueError(f'Unexpected A header: {path_a}')
    if list(b.columns) != ['ID','category','A','B','mass_times_Csum']:
        raise ValueError(f'Unexpected B header: {path_b}')
    for df, value in [(a, 'C_sum'), (b, 'mass_times_Csum')]:
        if df.isna().any().any() or (df.ID < 0).any():
            raise ValueError(f'Invalid/missing values in chunk {n}')
        if not df.A.between(0, 118).all() or not df.B.between(1, 511).all():
            raise ValueError(f'Invalid isotope codes in chunk {n}')
        if not np.isfinite(df[value]).all() or (df[value] < 0).any():
            raise ValueError(f'Invalid abundances in chunk {n}')
    ka, kb = encoded_keys(a), encoded_keys(b)
    ua, ub = np.unique(ka), np.unique(kb)
    if len(ua) != len(ka) or len(ub) != len(kb):
        raise ValueError(f'Duplicate ID/category/isotope keys in chunk {n}')
    if not np.array_equal(ua, ub):
        raise ValueError(f'A/B key mismatch in chunk {n}')
    same_order = np.array_equal(ka, kb)
    if not same_order:
        order_a, order_b = np.argsort(ka), np.argsort(kb)
        a = a.iloc[order_a].reset_index(drop=True)
        b = b.iloc[order_b].reset_index(drop=True)
    selected = a[a.category == args.category].copy()
    selected_b = b.loc[a.category == args.category].copy()
    species = selected.A.to_numpy(dtype=np.int64) * 512 + selected.B.to_numpy(dtype=np.int64)
    selected['species'] = species
    expected = sorted(key[0] * 512 + key[1] for col in columns
                      for kind, key in [request_kind(col)]
                      if kind == 'isotope' and col != 'X_{n}')
    if sorted(selected.species.unique().tolist()) != expected:
        raise ValueError(f'Unexpected/missing isotope set in chunk {n}')
    epochs = selected.groupby('ID')[['actual_timestep', 'actual_time']].nunique()
    if not (epochs == 1).all().all():
        raise ValueError(f'Inconsistent per-ID epoch in chunk {n}')
    category_counts = a.groupby(['ID', 'category'], observed=True).size()
    if len(category_counts) != 3 * a.ID.nunique() or not (category_counts == len(expected)).all():
        raise ValueError(f'Incomplete ID/category/species grid in chunk {n}')
    wide = selected.pivot(index='ID', columns='species', values='C_sum').sort_index()
    if wide.isna().any().any():
        raise ValueError(f'Incomplete selected isotope data in chunk {n}')
    ids = wide.index.to_numpy(dtype=np.int64)
    missing_hydro = ids[~np.isin(ids, physical.index)]
    hydro = physical.reindex(ids)
    peak = peaks.reindex(ids)
    output_values = {}
    for col in columns:
        kind, key = request_kind(col)
        if col == 'Id':
            output_values[col] = ids
        elif col in PHYSICAL:
            output_values[col] = hydro[PHYSICAL[col]].to_numpy()
        elif col in ['Peak temp', 'Peak density']:
            output_values[col] = peak[col].to_numpy()
        elif kind == 'isotope':
            code = key[0] * 512 + key[1]
            output_values[col] = wide[code].to_numpy() if code in wide else np.nan
        elif kind == 'element' and args.elements == 'subtotal':
            available = [c for c in wide if c // 512 == key[0]]
            output_values[col] = wide[available].sum(axis=1).to_numpy() if available else np.nan
        else:
            output_values[col] = np.nan
    output = pd.DataFrame(output_values, index=ids, columns=columns)
    # Check normalization for every ID/category, not a sum mixing three epochs.
    groups_a = a.groupby(['ID', 'category'], observed=True).C_sum.sum()
    groups_b = b.groupby(['ID', 'category'], observed=True).mass_times_Csum.sum()
    normalization = pd.DataFrame({'sum_C_sum': groups_a, 'sum_B_raw': groups_b}).reset_index()
    normalization['source_chunk'] = n
    normalization['AM_particle_mass_code'] = pd.to_numeric(
        physical.Mass.reindex(normalization.ID).reset_index(drop=True), errors='coerce')
    normalization['sum_B_solar_masses'] = normalization.sum_B_raw * 1e-6
    normalization['mass_fraction_residual'] = normalization.sum_C_sum - 1
    category_statistics = {}
    all_x = a.C_sum.to_numpy()
    all_b = b.mass_times_Csum.to_numpy()
    all_m = pd.to_numeric(physical.Mass.reindex(a.ID).reset_index(drop=True),
                          errors='coerce').to_numpy(dtype=np.float64)
    valid = np.isfinite(all_m)
    expected_b = all_m * all_x
    match_code = np.isclose(all_b, expected_b, rtol=2e-5, atol=1e-300)
    match_solar = np.isclose(all_b, expected_b * 1e-6, rtol=2e-5, atol=1e-300)
    for category in CATEGORIES:
        group = normalization[normalization.category == category]
        mask = (a.category == category).to_numpy()
        positive = mask & valid & (expected_b > 0)
        ratios = all_b[positive] / expected_b[positive]
        category_statistics[category] = {
            'IDs': len(group), 'sum_X_total': float(group.sum_C_sum.sum()),
            'sum_X_min': float(group.sum_C_sum.min()), 'sum_X_max': float(group.sum_C_sum.max()),
            'sum_X_eq_one_within_1e_5': int((group.mass_fraction_residual.abs() <= 1e-5).sum()),
            'sum_X_eq_one_within_1e_4': int((group.mass_fraction_residual.abs() <= 1e-4).sum()),
            'sum_X_eq_one_within_1e_3': int((group.mass_fraction_residual.abs() <= 1e-3).sum()),
            'sum_X_zero': int((group.sum_C_sum == 0).sum()),
            'sum_X_above_1_plus_1e_4': int((group.sum_C_sum > 1.0001).sum()),
            'sum_X_below_1_minus_1e_4': int((group.sum_C_sum < 0.9999).sum()),
            'sum_B_raw': float(group.sum_B_raw.sum()),
            'sum_B_solar_masses': float(group.sum_B_solar_masses.sum()),
            'sum_AM_particle_mass_code': float(group.AM_particle_mass_code.sum()),
            'sum_AM_particle_mass_solar': float(group.AM_particle_mass_code.sum() * 1e-6),
            'B_code_unit_comparisons': int((mask & valid).sum()),
            'B_code_unit_mismatches': int((mask & valid & ~match_code).sum()),
            'B_positive_comparisons': int(positive.sum()),
            'B_positive_matches_if_already_solar': int((positive & match_solar).sum()),
            'B_divided_by_code_mass_times_X_min': float(ratios.min()) if len(ratios) else None,
            'B_divided_by_code_mass_times_X_max': float(ratios.max()) if len(ratios) else None}
    # The B values are a separate mass-weighted quantity, never substitutes for X.
    x = selected.C_sum.to_numpy()
    weighted = selected_b.mass_times_Csum.to_numpy()
    masses = pd.to_numeric(physical.Mass.reindex(selected.ID).reset_index(drop=True),
                           errors='coerce').to_numpy(dtype=np.float64)
    comparable = np.isfinite(masses)
    expected_weight = masses * x
    mass_match = np.isclose(weighted, expected_weight, rtol=2e-5, atol=1e-300)
    relative = np.zeros(len(x))
    nonzero = comparable & (expected_weight != 0)
    relative[nonzero] = np.abs(weighted[nonzero] / expected_weight[nonzero] - 1)
    mismatch = comparable & ~mass_match
    sample = []
    for j in np.flatnonzero(mismatch)[:5]:
        row = selected.iloc[j]
        sample.append({'Id': int(row.ID), 'Z': int(row.A), 'A': int(row.B),
                       'C_sum': float(x[j]), 'B_value': float(weighted[j]),
                       'AM_mass_times_C_sum': float(expected_weight[j])})
    epoch = selected.drop_duplicates('ID').set_index('ID').reindex(ids)
    sidecar = pd.DataFrame({'Id': ids, 'source_chunk': n, 'category': args.category,
                           'actual_timestep': epoch.actual_timestep.to_numpy(),
                           'actual_time_as_written': epoch.actual_time.to_numpy(),
                           'hydro_snapshot_iteration': 2000,
                           'hydro_snapshot_count': peak.hydro_snapshot_count.to_numpy(),
                           'hydro_matched': np.isin(ids, physical.index),
                           'selected_isotope_sum': wide.sum(axis=1).to_numpy(),
                           'source_zero_isotope_count': (wide == 0).sum(axis=1).to_numpy()})
    part = WORK / 'parts' / f'{n:04d}.out'
    output.to_csv(part, sep='\t', index=False, na_rep='NA', float_format='%.10g', mode='x')
    sidecar.to_csv(WORK / 'parts' / f'{n:04d}.epochs.csv', index=False,
                   na_rep='NA', float_format='%.17g', mode='x')
    normalization.to_csv(WORK / 'parts' / f'{n:04d}.normalization.csv', index=False,
                         na_rep='NA', float_format='%.17g', mode='x')
    record = {'chunk': n, 'rows_A': len(a), 'rows_B': len(b),
              'category': args.category, 'element_policy': args.elements,
              'category_statistics': category_statistics,
              'IDs': len(ids), 'A_B_same_row_order': same_order,
              'A_B_keys_match': True, 'duplicate_keys': 0,
              'missing_physical_IDs': missing_hydro.tolist(),
              'selected_B_mass_comparisons': int(comparable.sum()),
              'selected_B_mass_mismatches': int(mismatch.sum()),
              'selected_B_mass_max_relative_error': float(relative[nonzero].max()) if nonzero.any() else 0.,
              'selected_B_mass_mismatch_examples': sample,
              'selected_isotope_sum_min': float(wide.sum(axis=1).min()),
              'selected_isotope_sum_max': float(wide.sum(axis=1).max()),
              'selected_time_min': float(epoch.actual_time.min()),
              'selected_time_max': float(epoch.actual_time.max()),
              'source_A': meta_a, 'source_B': meta_b}
    save_json(WORK / 'parts' / f'{n:04d}.json', record)
    return record, ids


def assemble(args, columns, records, physical_check, peaks_meta, all_ids):
    output_path = ROOT / 'final_chunks.out'
    temp = WORK / 'final_chunks.out.partial'
    times_temp = WORK / 'grain_epochs.csv.partial'
    norm_temp = WORK / 'normalization_by_id.csv.partial'
    missing_counts = dict.fromkeys(columns, 0)
    seen = set()
    digest = hashlib.sha256()
    with temp.open('xb') as out, times_temp.open('xb') as times, norm_temp.open('xb') as norms:
        for n in sorted(records):
            part = WORK / 'parts' / f'{n:04d}.out'
            with part.open('rb') as inp:
                header = inp.readline()
                if header.decode().rstrip('\n').split('\t') != columns:
                    raise ValueError(f'Incorrect output columns in chunk {n}')
                if n == min(records):
                    out.write(header); digest.update(header)
                count = 0
                for line in inp:
                    values = line.rstrip(b'\n').split(b'\t')
                    if len(values) != len(columns):
                        raise ValueError(f'Incorrect output row width in chunk {n}')
                    ident = int(values[0])
                    if ident in seen:
                        raise ValueError(f'Duplicate output ID: {ident}')
                    seen.add(ident)
                    for col, value in zip(columns, values):
                        if value == b'NA': missing_counts[col] += 1
                        elif not value: raise ValueError(f'Unmarked empty value at ID {ident}/{col}')
                    out.write(line); digest.update(line); count += 1
            if count != records[n]['IDs']:
                raise ValueError(f'Partition row-count mismatch: chunk {n}')
            with (WORK / 'parts' / f'{n:04d}.epochs.csv').open('rb') as inp:
                h = inp.readline()
                if n == min(records): times.write(h)
                while True:
                    block = inp.read(1024 * 1024)
                    if not block: break
                    times.write(block)
            with (WORK / 'parts' / f'{n:04d}.normalization.csv').open('rb') as inp:
                h = inp.readline()
                if n == min(records): norms.write(h)
                while True:
                    block = inp.read(1024 * 1024)
                    if not block: break
                    norms.write(block)
            if n % 100 == 0: log(f'Output verification: {n}/{len(records)} chunks')
    if seen != set(all_ids):
        raise ValueError('Final ID set differs from validated source IDs')
    if output_path.exists():
        raise FileExistsError(output_path)
    # Hard link publishes the validated artifact atomically and refuses overwrite.
    os.link(temp, output_path)
    temp.unlink()  # Only our own temporary output, never an input.
    os.link(times_temp, WORK / 'grain_epochs.csv')
    times_temp.unlink()
    os.link(norm_temp, WORK / 'normalization_by_id.csv')
    norm_temp.unlink()
    category_totals = {}
    for category in CATEGORIES:
        stats = [r['category_statistics'][category] for r in records.values()]
        totals = {}
        for key in stats[0]:
            values = [s[key] for s in stats if s[key] is not None]
            if key.endswith('_min'): totals[key] = min(values) if values else None
            elif key.endswith('_max'): totals[key] = max(values) if values else None
            else: totals[key] = sum(values)
        category_totals[category] = totals
    mass_report = {'category_totals': category_totals,
                   'B_raw_units': 'native code mass; 1 code mass = 1e-6 solar masses under the local SNSPH convention',
                   'unit_evidence': [
                       '/data/pyoung3/gsvance/data_mining/sn_utils.py:58: SNSPH_MASS = 1e-6 * 1.9889e33',
                       '/data/pyoung3/snsphtools/SDFtoASCII.c:156: massCF default 1.9889e27 grams',
                       '/scratch/svierlin/prism-v02/sum_eoc_manifest_part2.sh: mass * C_sum without conversion',
                       '/scratch/svierlin/prism-v02/sum_eoc_manifest_part3_RP.sh:251: sum[k] * 1e-6'],
                   'unit_provenance': 'local source convention; AM SDF headers omit massCF',
                   'interpretation': 'Per-category B sums cover the supplied IDs and selected isotopes only; they do not independently establish the total ejected mass. Do not add the three categories together.'}
    save_json(WORK / 'mass_fraction_validation.json', mass_report)
    pd.DataFrame.from_dict(category_totals, orient='index').rename_axis('category').to_csv(
        WORK / 'mass_fraction_summary.csv', float_format='%.17g', mode='x')
    summary = {'output': str(output_path), 'rows': len(seen), 'columns': len(columns),
               'header_exactly_matches_requested_list': True, 'delimiter': 'tab',
               'output_bytes': output_path.stat().st_size, 'output_sha256': digest.hexdigest(),
               'category': args.category, 'element_policy': args.elements,
               'chunk_pairs': len(records), 'duplicate_output_IDs': 0,
               'missing_cells_by_column': missing_counts,
               'physical_verification': physical_check,
               'peak_snapshot_count': len(peaks_meta['snapshots']),
               'B_mass_mismatches': sum(r['selected_B_mass_mismatches'] for r in records.values()),
               'B_mass_comparisons': sum(r['selected_B_mass_comparisons'] for r in records.values()),
               'B_mass_relative_tolerance': 2e-5,
               'category_totals': category_totals,
               'source_rows_A': sum(r['rows_A'] for r in records.values()),
               'source_rows_B': sum(r['rows_B'] for r in records.values()),
               'actual_time_min': min(r['selected_time_min'] for r in records.values()),
               'actual_time_max': max(r['selected_time_max'] for r in records.values()),
               'physical_IDs_without_chunks': len(set(physical_check.pop('_physical_ids', [])) - seen),
               'source_metadata_unchanged_at_finish': True}
    summary['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    summary['requested_columns_sha256'] = hashlib.sha256(
        (ROOT / 'List of isotopes_elements_variables.txt').read_bytes()).hexdigest()
    summary['numpy_version'] = np.__version__
    summary['pandas_version'] = pd.__version__
    save_json(WORK / 'validation.json', summary)
    save_json(WORK / 'chunk_validation.json', [records[n] for n in sorted(records)])
    unavailable = [name for name, count in missing_counts.items() if count == len(seen)]
    notes = f'''AM chunk export: final_chunks.out

Rows: {len(seen):,}; columns: {len(columns)}; delimiter: TAB; missing value: NA.
One header row, exactly the order/spelling in List of isotopes_elements_variables.txt.
One row per ID represented in the supplied chunk files; no invented IDs.
Rows are ordered by numeric chunk number, then numeric ID within each chunk.
Id is the source simulation/SPH particle ID (the user's grain/clump identifier).

Composition category: {args.category}. Isotope columns copy C_sum as dimensionless
mass-fraction values from final_chunks/chunk_N.csv, without normalization or decay.
The two chunk sets were joined/checked using ID, category, proton number (source A),
and mass number (source B). B contains particle-mass-weighted isotope values;
these are checked separately and are NOT substituted for isotope mass fractions.
Mass comparison tolerance is 2e-5 relative because source products are rounded.
B mass discrepancies: {summary['B_mass_mismatches']} of {summary['B_mass_comparisons']} checked cells.
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
Y_{{e}} retains the stored Y_e value; this is not a recomputed/validated electron
fraction. No CGS conversion is applied. Hydro iteration is 2000 and tpos is
1575.58 code time. Composition actual_time is retained separately in grain_epochs.csv;
matching IDs does not imply matching physical epochs. Historical trajectory time
and density problems described in the existing audit are not repaired by reshaping.

Peak temp is the maximum stored temp (K) and Peak density the maximum stored rho
(code density) across {len(peaks_meta['snapshots'])} saved raw AM snapshots. These are separate
maxima, not necessarily simultaneous, and are sampled hydro maxima, not guaranteed
continuous peaks or maxima of the PRISM input trajectories. Per-ID snapshot counts
are in grain_epochs.csv; complete source coverage is in peak_sources.json.

Element policy: {args.elements}. With 'na', complete elemental totals are unknown
because chunks contain only selected isotopes. With 'subtotal', element columns
sum only the available selected isotopes and must be interpreted as subtotals.
An element with no selected isotope remains NA. X_{{n}} is absent from the chunks
and remains NA; hydro-network neutron fractions are not mixed into PRISM composition.
Columns entirely unavailable: {', '.join(unavailable)}.

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
'''
    (WORK / 'README.txt').write_text(notes)
    log(f'COMPLETE: {output_path}; {len(seen):,} rows, {len(columns)} columns, sha256={digest.hexdigest()}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--category', choices=CATEGORIES, default='last')
    parser.add_argument('--elements', choices=['na', 'subtotal'], default='na')
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--limit', type=int, default=938)
    args = parser.parse_args()
    (WORK / 'parts').mkdir(exist_ok=True)
    columns = (ROOT / 'List of isotopes_elements_variables.txt').read_text().splitlines()
    if len(columns) != len(set(columns)) or columns[0] != 'Id':
        raise ValueError('Invalid requested column list')
    output_path = ROOT / 'final_chunks.out'
    if output_path.exists(): raise FileExistsError(output_path)
    for folder, pattern, regex in [
        (ROOT / 'final_chunks', 'chunk_*.csv', r'chunk_(\d+)\.csv'),
        (ROOT / 'final_chunks/final_chunks_B', 'chunk_*_B.csv', r'chunk_(\d+)_B\.csv')]:
        numbers = sorted(int(re.fullmatch(regex, p.name)[1]) for p in folder.glob(pattern))
        if numbers != list(range(1, 939)):
            raise ValueError(f'Unexpected/missing chunk numbers in {folder}')
    physical, physical_check = get_physical()
    peaks, peaks_meta = get_peaks()
    records, all_ids = {}, []
    log(f'Processing {args.limit} chunk pairs with {args.workers} workers')
    started = time.monotonic()
    with futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending = {}
        for n in range(1, args.limit + 1):
            done_path = WORK / 'parts' / f'{n:04d}.json'
            if done_path.exists():
                record = json.loads(done_path.read_text())
                if record.get('category') != args.category or record.get('element_policy') != args.elements:
                    raise ValueError(f'Cached chunk {n} has different export options')
                for key in ['source_A', 'source_B']:
                    if fingerprint(Path(record[key]['path'])) != {k:record[key][k] for k in ['path','bytes','mtime_ns']}:
                        raise ValueError(f'Changed source for cached chunk {n}')
                ids = pd.read_csv(WORK / 'parts' / f'{n:04d}.epochs.csv', usecols=['Id']).Id.to_numpy()
                records[n] = record; all_ids.extend(ids.tolist())
            else:
                pending[pool.submit(process_chunk, n, args, columns, physical, peaks)] = n
        for future in futures.as_completed(pending):
            try:
                record, ids = future.result()
            except Exception:
                for other in pending: other.cancel()
                raise
            records[record['chunk']] = record
            all_ids.extend(ids.tolist())
            if len(records) % 20 == 0 or len(records) == args.limit:
                log(f'Chunks: {len(records)}/{args.limit}; {len(all_ids):,} IDs; elapsed {time.monotonic()-started:.0f}s')
    if len(set(all_ids)) != len(all_ids): raise ValueError('IDs duplicated across chunks')
    # Verify source stat fingerprints once again before publishing any combined file.
    sources = [r[key] for r in records.values() for key in ['source_A', 'source_B']]
    sources += [physical_check['source'], physical_check['raw_snapshot']]
    sources += [r['source'] for r in peaks_meta['snapshots']]
    for record in sources:
        if fingerprint(Path(record['path'])) != {k:record[k] for k in ['path','bytes','mtime_ns']}:
            raise ValueError(f'Source changed during export: {record["path"]}')
    if args.limit != 938:
        log('Pilot complete; final combined output is intentionally not published')
        return
    physical_check['_physical_ids'] = physical.index.tolist()
    missing = physical.loc[~physical.index.isin(all_ids), ['ID', 'Mass']].copy()
    missing.to_csv(WORK / 'hydro_ids_without_chunks.csv', index=False, mode='x')
    physical_check['all_hydro_mass_code'] = float(pd.to_numeric(physical.Mass).sum())
    physical_check['hydro_mass_without_chunks_code'] = float(pd.to_numeric(missing.Mass).sum())
    physical_check['hydro_IDs_without_chunks'] = len(missing)
    assemble(args, columns, records, physical_check, peaks_meta, all_ids)


if __name__ == '__main__':
    main()
