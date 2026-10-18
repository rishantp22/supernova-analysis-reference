#!/usr/bin/env python3
"""Read-only, bounded-memory audit of two isotope tables; outputs in a new directory.

Uses disk-backed GNU sort and streaming NumPy blocks, never a whole-file DataFrame.
All matching uses integer particle IDs and normalized column names.
"""
import argparse
from collections import Counter
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import resource
import subprocess
import time
import warnings

import numpy as np

FACTORS = {**dict.fromkeys(['x', 'y', 'z', 'h'], 6.955e10),
           **dict.fromkeys(['vx', 'vy', 'vz'], 6.955e8),
           **dict.fromkeys(['ax', 'ay', 'az'], 6.955e6),
           'mass': 1.9889e27,
           'density': 1.9889e27 / 6.955e10**3,
           'peak density': 1.9889e27 / 6.955e10**3}
MISSING = b'NA'
BLOCK = 1024


def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True)


def fingerprint(path):
    s = path.stat()
    return dict(path=str(path.resolve()), bytes=s.st_size, mtime_ns=s.st_mtime_ns,
                ctime_ns=s.st_ctime_ns, inode=s.st_ino, device=s.st_dev)


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def canonical(name):
    if name == 'Y_{e}':
        return name
    if name.startswith('X_{') and name.endswith('}'):
        return name[3:-1]
    return name.lower() if name.lower() in ['id', 'density', 'peak temp', 'peak density'] else name


def numbers(rows, ncols):
    payload = b'\n'.join(rows)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        result = np.fromstring(payload.replace(b'NA', b'nan'), dtype=np.float64, sep=' ')
    if result.size != len(rows) * ncols:
        raise ValueError('Invalid numeric token or unexpected numeric column count')
    result = result.reshape(len(rows), ncols)
    if np.isinf(result).any() or int(np.isnan(result).sum()) != sum(r.split().count(MISSING) for r in rows):
        raise ValueError('Unexpected NaN/infinity or malformed missing-value token')
    return result


def audited_numbers(rows, ncols):
    """Keep invalid numeric tokens distinct from the explicit NA marker."""
    try:
        result = numbers(rows, ncols)
        return result, np.zeros(result.shape, dtype=bool)
    except (ValueError, Warning):
        result = np.full((len(rows), ncols), np.nan)
        invalid = np.zeros(result.shape, dtype=bool)
        for i, row in enumerate(rows):
            tokens = row.split(b'\t')
            if len(tokens) != ncols:
                raise ValueError('Unexpected width in already validated row')
            for j, token in enumerate(tokens):
                if token == MISSING:
                    continue
                try:
                    value = float(token)
                    if not math.isfinite(value):
                        raise ValueError()
                    result[i, j] = value
                except ValueError:
                    invalid[i, j] = True
        return result, invalid


def scan(path, delimiter, work, tag):
    log(f'Scanning {path.name}')
    before = fingerprint(path)
    digest = hashlib.sha256()
    nonblank_digest = hashlib.sha256()
    sequence = hashlib.sha256()
    counts = Counter()
    malformed = []
    blank = invalid_ids = rows = inversions = equal_adjacent = 0
    first_ids, last_ids = [], []
    previous = None
    payloads = []
    payload_lines = []
    invalid_examples, invalid_id_examples = [], []
    keyed = work / f'{tag}.keyed.tsv'
    with path.open('rb', buffering=4*1024*1024) as source, keyed.open('xb', buffering=4*1024*1024) as output:
        header_raw = source.readline()
        digest.update(header_raw)
        nonblank_digest.update(header_raw)
        header = [s.strip().decode() for s in header_raw.rstrip(b'\r\n').split(delimiter)]
        names = [canonical(s) for s in header]
        if names[0] != 'id' or len(names) != len(set(names)):
            raise ValueError('Ambiguous normalized header or missing first ID column')
        n = len(names)
        missing = np.zeros(n, np.int64)
        zeros = np.zeros(n, np.int64)
        negative = np.zeros(n, np.int64)
        invalid_counts = np.zeros(n, np.int64)
        numeric_rows = 0
        invalid_numeric_rows = 0

        def flush():
            nonlocal numeric_rows, invalid_numeric_rows
            if not payloads:
                return
            values, invalid = audited_numbers(payloads, n)
            missing[:] += (np.isnan(values) & ~invalid).sum(axis=0)
            invalid_counts[:] += invalid.sum(axis=0)
            invalid_numeric_rows += int(invalid.any(axis=1).sum())
            for i,j in zip(*np.nonzero(invalid)):
                if len(invalid_examples) < 20:
                    tokens = payloads[i].split(b'\t')
                    invalid_examples.append(dict(line=payload_lines[i], id=tokens[0].decode(),
                                                 column=header[j], token=tokens[j].decode(errors='replace')))
            zeros[:] += (values == 0).sum(axis=0)
            negative[:] += (values < 0).sum(axis=0)
            numeric_rows += len(payloads)
            payloads.clear()
            payload_lines.clear()

        last_raw = header_raw
        for line_number, raw in enumerate(source, 2):
            digest.update(raw)
            last_raw = raw
            if not raw.strip():
                blank += 1
                continue
            nonblank_digest.update(raw)
            rows += 1
            # Both observed formats are unquoted numeric tables with explicit NA.
            tokens = [s.strip() for s in raw.rstrip(b'\r\n').split(delimiter)]
            counts[len(tokens)] += 1
            if len(tokens) != n:
                if len(malformed) < 10:
                    malformed.append(dict(line=line_number, fields=len(tokens), id=tokens[0].decode(errors='replace')))
                continue
            if not re.fullmatch(rb'[0-9]+', tokens[0]):
                invalid_ids += 1
                if len(invalid_id_examples) < 20:
                    invalid_id_examples.append(dict(line=line_number, token=tokens[0].decode(errors='replace')))
                continue
            ident = int(tokens[0])
            if previous is not None:
                inversions += ident < previous
                equal_adjacent += ident == previous
            previous = ident
            if len(first_ids) < 10:
                first_ids.append(ident)
            last_ids.append(ident)
            last_ids[:] = last_ids[-10:]
            sequence.update(str(ident).encode() + b'\n')
            normalized = b'\t'.join(tokens)
            output.write(str(ident).encode() + b'\t' + str(line_number).encode() + b'\t' + b'\t'.join(tokens[1:]) + b'\n')
            payloads.append(normalized)
            payload_lines.append(line_number)
            if len(payloads) == BLOCK:
                flush()
            if rows % 100000 == 0:
                log(f'{tag}: scanned {rows:,} rows')
        flush()
    after = fingerprint(path)
    if before != after:
        raise RuntimeError(f'Input changed while being read: {path}')
    stats = dict(metadata=before, sha256=digest.hexdigest(), nonblank_sha256=nonblank_digest.hexdigest(), header=header,
                 normalized_header=names, columns=n, rows=rows, numeric_rows=numeric_rows,
                 field_counts=dict(counts), malformed_examples=malformed,
                 malformed_rows=rows-numeric_rows-invalid_ids, blank_lines=blank,
                 invalid_or_missing_ids=invalid_ids, numeric_order_inversions=inversions,
                 invalid_id_examples=invalid_id_examples,
                 invalid_numeric_cells=int(invalid_counts.sum()), invalid_numeric_rows=invalid_numeric_rows,
                 invalid_numeric_by_column=dict(zip(names,invalid_counts.tolist())),
                 invalid_numeric_examples=invalid_examples,
                 adjacent_equal_ids=equal_adjacent, first_ids=first_ids, last_ids=last_ids,
                 id_sequence_sha256=sequence.hexdigest(), ends_with_newline=last_raw.endswith(b'\n'),
                 missing_by_column=dict(zip(names, missing.tolist())),
                 zeros_by_column=dict(zip(names, zeros.tolist())),
                 negative_by_column=dict(zip(names, negative.tolist())), input_unchanged_during_scan=True)
    save(work / f'{tag}_scan.json', stats)
    if inversions:
        ordered = work / f'{tag}.sorted.tsv'
        log(f'{tag}: sorting temporary records by numeric particle ID (128 MiB sort buffer)')
        subprocess.run(['sort', '-t', '\t', '-k1,1n', '-s', '--buffer-size=128M', '--parallel=1',
                        '-T', str(work), '-o', str(ordered), str(keyed)], check=True,
                       env={**os.environ, 'LC_ALL': 'C'})
        keyed.unlink()  # Only our derived temporary file; never an input.
    else:
        ordered = keyed
    return stats, ordered


def grouped(path, stats, duplicate_writer):
    unique = dup_ids = extra = 0
    with path.open('rb', buffering=4*1024*1024) as stream:
        for key, group in itertools.groupby(stream, key=lambda line: int(line.split(b'\t', 1)[0])):
            first = next(group)
            count = 1
            sample_lines = [int(first.split(b'\t', 2)[1])]
            for line in group:
                count += 1
                if len(sample_lines) < 10:
                    sample_lines.append(int(line.split(b'\t', 2)[1]))
            unique += 1
            if count > 1:
                dup_ids += 1
                extra += count-1
                duplicate_writer.writerow([key, count, ' '.join(map(str, sample_lines))])
            yield key, count, first.rstrip(b'\n').split(b'\t', 2)[2]
    stats.update(unique_ids=unique, duplicate_ids=dup_ids, duplicate_extra_rows=extra)


def compare(a, b, apath, bpath, work):
    log('Matching sorted records by ID and columns by normalized names')
    names = a['normalized_header'][1:]
    bnames = b['normalized_header'][1:]
    if set(names) != set(bnames):
        raise ValueError('Column sets differ; inspect saved headers before comparison')
    reorder = [bnames.index(name) for name in names]
    factors = np.array([FACTORS.get(name, 1.) for name in names])
    isotope = np.array([bool(re.fullmatch(r'[0-9]+[A-Z][a-z]?', name)) for name in names])
    summary = dict(shared_unique_ids=0, only_final_chunks_ids=0, only_plotting_ids=0,
                   ambiguous_shared_duplicate_ids=0, compared_unique_ids=0,
                   isotope_columns=int(isotope.sum()), composition_text_identical_rows=0,
                   rows_with_raw_numeric_differences=0, rows_with_normalized_differences_1e9=0,
                   rows_with_isotope_differences=0)
    metrics = {key: np.zeros(len(names), np.int64) for key in
               ['finite_pairs', 'both_missing', 'missing_mismatches', 'invalid_numeric_pairs', 'raw_exact_differences',
                'converted_exact_differences', 'converted_mismatches_rtol_1e9',
                'converted_mismatches_rtol_1e5']}
    maxabs = np.zeros(len(names))
    maxrel = np.zeros(len(names))
    examples = {name: [] for name in names}
    blocks_a, blocks_b, ids = [], [], []
    sample = {}

    def flush():
        if not ids:
            return
        av, ia = audited_numbers(blocks_a, len(names))
        bv, ib = audited_numbers(blocks_b, len(bnames))
        bv, ib = bv[:,reorder], ib[:,reorder]
        invalid = ia | ib
        ma, mb = np.isnan(av) & ~ia, np.isnan(bv) & ~ib
        finite = ~ma & ~mb & ~invalid
        mismatch_na = ma != mb
        converted = av * factors
        absolute = np.zeros_like(av)
        np.subtract(converted, bv, out=absolute, where=finite)
        np.abs(absolute, out=absolute)
        relative = np.zeros_like(av)
        np.divide(absolute, np.maximum(np.abs(converted), np.abs(bv)), out=relative,
                  where=finite & ((converted != 0) | (bv != 0)))
        exact_raw = finite & (av != bv)
        exact_converted = finite & (converted != bv)
        # Explicitly use B as the reference in this documented asymmetric formula.
        strict = finite & (absolute > (1e-300 + 1e-9 * np.abs(bv)))
        printed = finite & (absolute > (1e-300 + 1e-5 * np.abs(bv)))
        for key, mask in [('finite_pairs', finite), ('both_missing', ma & mb),
                          ('invalid_numeric_pairs', invalid),
                          ('missing_mismatches', mismatch_na), ('raw_exact_differences', exact_raw),
                          ('converted_exact_differences', exact_converted),
                          ('converted_mismatches_rtol_1e9', strict),
                          ('converted_mismatches_rtol_1e5', printed)]:
            metrics[key] += mask.sum(axis=0)
        maxabs[:] = np.maximum(maxabs, absolute.max(axis=0))
        maxrel[:] = np.maximum(maxrel, relative.max(axis=0))
        summary['rows_with_raw_numeric_differences'] += int((exact_raw | mismatch_na | invalid).any(axis=1).sum())
        summary['rows_with_normalized_differences_1e9'] += int((strict | mismatch_na | invalid).any(axis=1).sum())
        summary['rows_with_isotope_differences'] += int((exact_raw[:, isotope] | mismatch_na[:, isotope] | invalid[:,isotope]).any(axis=1).sum())
        for j in np.flatnonzero((strict | mismatch_na).any(axis=0)):
            name = names[j]
            for i in np.flatnonzero(strict[:, j] | mismatch_na[:, j])[:max(0, 5-len(examples[name]))]:
                examples[name].append(dict(id=ids[i], final_chunks=None if ma[i,j] else float(av[i,j]),
                                           plotting=None if mb[i,j] else float(bv[i,j]),
                                           converted_final_chunks=None if ma[i,j] else float(converted[i,j]),
                                           absolute_difference=float(absolute[i,j]),
                                           symmetric_relative_difference=float(relative[i,j])))
        for particle in [23, 7753, 10000]:
            if particle in ids:
                i = ids.index(particle)
                sample[str(particle)] = {name: dict(final_chunks=float(av[i,j]), plotting=float(bv[i,j]),
                                                    converted_final_chunks=float(converted[i,j]))
                                        for j, name in enumerate(names) if j < 15 or name in ['4He','12C','16O','56Fe','56Ni']}
        summary['compared_unique_ids'] += len(ids)
        if summary['compared_unique_ids'] % (BLOCK*100) == 0:
            log(f"Compared {summary['compared_unique_ids']:,} matched IDs")
        blocks_a.clear(); blocks_b.clear(); ids.clear()

    with (work/'duplicate_ids_final_chunks.csv').open('x') as da, (work/'duplicate_ids_plotting.csv').open('x') as db, \
         (work/'ids_only_final_chunks.csv').open('x') as oa, (work/'ids_only_plotting.csv').open('x') as ob:
        wa, wb = csv.writer(da), csv.writer(db)
        wa.writerow(['id','occurrences','first_source_lines']); wb.writerow(['id','occurrences','first_source_lines'])
        xa, xb = csv.writer(oa), csv.writer(ob)
        xa.writerow(['id']); xb.writerow(['id'])
        ga, gb = grouped(apath, a, wa), grouped(bpath, b, wb)
        ra, rb = next(ga, None), next(gb, None)
        while ra is not None or rb is not None:
            if rb is None or (ra is not None and ra[0] < rb[0]):
                summary['only_final_chunks_ids'] += 1; xa.writerow([ra[0]])
                ra = next(ga, None)
            elif ra is None or rb[0] < ra[0]:
                summary['only_plotting_ids'] += 1; xb.writerow([rb[0]])
                rb = next(gb, None)
            else:
                summary['shared_unique_ids'] += 1
                if ra[1] != 1 or rb[1] != 1:
                    summary['ambiguous_shared_duplicate_ids'] += 1
                else:
                    ids.append(ra[0]); blocks_a.append(ra[2]); blocks_b.append(rb[2])
                    if names == bnames and ra[2].split(b'\t', 15)[15] == rb[2].split(b'\t', 15)[15]:
                        summary['composition_text_identical_rows'] += 1
                    if len(ids) == BLOCK:
                        flush()
                ra, rb = next(ga, None), next(gb, None)
        flush()
    columns = []
    for j, name in enumerate(names):
        columns.append(dict(column=name, final_chunks_name=a['header'][j+1],
                            plotting_name=b['header'][reorder[j]+1], conversion_factor=float(factors[j]),
                            is_isotope=bool(isotope[j]), **{k:int(v[j]) for k,v in metrics.items()},
                            max_converted_absolute_difference=float(maxabs[j]),
                            max_symmetric_relative_difference=float(maxrel[j])))
    with (work/'column_comparison.csv').open('x') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(columns[0])); writer.writeheader(); writer.writerows(columns)
    summary['isotope_finite_pairs'] = int(metrics['finite_pairs'][isotope].sum())
    summary['isotope_exact_numeric_mismatches'] = int(metrics['raw_exact_differences'][isotope].sum())
    summary['isotope_missing_mismatches'] = int(metrics['missing_mismatches'][isotope].sum())
    summary['all_column_missing_mismatches'] = int(metrics['missing_mismatches'].sum())
    summary['invalid_numeric_pairs'] = int(metrics['invalid_numeric_pairs'].sum())
    summary['column_order_equal_after_name_mapping'] = names == bnames
    summary['row_order_identical'] = a['id_sequence_sha256'] == b['id_sequence_sha256']
    save(work/'mismatch_examples.json', {k:v for k,v in examples.items() if v})
    save(work/'selected_particle_examples.json', sample)
    return summary, columns


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--final', type=Path, required=True)
    parser.add_argument('--plotting', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-final-bytes', type=int, default=0,
                        help='Minimum complete length from prior manifest; differences are audited, not rejected')
    parser.add_argument('--expected-final-sha256')
    parser.add_argument('--expected-plotting-sha256')
    parser.add_argument('--wait-minutes', type=float, default=30)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    deadline = start + args.wait_minutes*60
    previous, stable = None, 0
    while True:
        now = [fingerprint(args.final), fingerprint(args.plotting)]
        ready = not args.expected_final_bytes or now[0]['bytes'] >= args.expected_final_bytes
        stable = stable+1 if ready and now == previous else 0
        if stable >= 3:
            break
        if time.monotonic() > deadline:
            raise RuntimeError('Input did not reach expected complete size and remain stable for 30 seconds')
        log(f"Waiting for stable complete input: {now[0]['bytes']:,} / {args.expected_final_bytes:,} bytes")
        previous = now
        time.sleep(10)
    before = now
    byte_result = subprocess.run(['cmp', str(args.final), str(args.plotting)], capture_output=True, text=True)
    if byte_result.returncode not in (0, 1):
        raise RuntimeError(byte_result.stderr)
    a, apath = scan(args.final, b'\t', args.output, 'final_chunks')
    b, bpath = scan(args.plotting, b',', args.output, 'plotting')
    summary, columns = compare(a, b, apath, bpath, args.output)
    after = [fingerprint(args.final), fingerprint(args.plotting)]
    if before != after:
        raise RuntimeError('An original changed during the audit; results invalid')
    for info, expected in [(a,args.expected_final_sha256),(b,args.expected_plotting_sha256)]:
        info['matches_prior_export_sha256'] = info['sha256'] == expected if expected else None
        info['matches_prior_export_after_ignoring_blank_lines'] = info['nonblank_sha256'] == expected if expected else None
    report = dict(byte_identical=byte_result.returncode == 0, cmp_message=byte_result.stdout.strip(),
                  final_chunks=a, plotting=b, comparison=summary, columns=columns,
                  tolerance=dict(exact_isotope_comparison='float64 equality plus normalized text equality',
                                 formula='abs(final * factor - plotting) <= atol + rtol * abs(plotting)',
                                 strict_rtol=1e-9, printed_precision_rtol=1e-5, atol=1e-300,
                                 missing_policy='NA is missing, never zero; matching NA counted separately'),
                  source_metadata_unchanged=True, elapsed_seconds=time.monotonic()-start,
                  max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                  unit_provenance='/scratch/rpraka15/final_chunks_export_20261018/make_plotting_format.py')
    save(args.output/'comparison.json', report)
    log(json.dumps(summary))
    log(f"COMPLETE; peak Python RSS {report['max_rss_kib']/1024:.1f} MiB; originals unchanged")
    apath.unlink(); bpath.unlink()  # Remove only generated sort intermediates.


if __name__ == '__main__':
    main()
