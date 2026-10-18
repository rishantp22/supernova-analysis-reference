#!/usr/bin/env python3
"""Names-only particle census; never stat particle entries or recurse into them."""
import argparse
import csv
import datetime
import json
import os
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('/data/pyoung3/groupdata/svierlin'))
    parser.add_argument('--output', type=Path, default=Path('/scratch/rpraka15/supernova-grain-analysis/inventory'))
    args = parser.parse_args()
    args.output = args.output.resolve()
    allowed = [Path('/scratch/rpraka15/supernova-grain-analysis').resolve(),
               Path('/home/rpraka15/supernova-grain-analysis').resolve()]
    if not any(args.output == root or root in args.output.parents for root in allowed):
        parser.error('--output must be inside the user-owned supernova-grain-analysis home or scratch project')
    args.output.mkdir(parents=True, exist_ok=True)
    report = {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scope': 'Direct entry names only; numeric names are candidate particle IDs, not verified directories or completed runs.',
              'models': {}, 'archives': [], 'errors': []}
    samples = {}
    for model, subdir in [('jet3b', 'jet3b_data/jet3b_outputs'), ('cco2', 'cco2_data/cco2_outputs'), ('vconv', 'vconv_data/vconv_outputs')]:
        source = args.source / subdir
        ids = []
        other = []
        try:
            with os.scandir(source) as entries:
                for entry in entries:
                    if entry.name.isascii() and entry.name.isdecimal():
                        ids.append(entry.name)
                    else:
                        other.append(entry.name)
        except OSError as exc:
            report['errors'].append({'source': str(source), 'error': str(exc)})
        ids.sort(key=lambda value: (int(value), value))
        with (args.output / f'{model}_particle_ids.csv').open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['model', 'particle_id', 'source_path'])
            writer.writerows((model, pid, str(source / pid)) for pid in ids)
        k = min(24, len(ids))
        indices = [round(i * (len(ids) - 1) / (k - 1)) for i in range(k)] if k > 1 else ([0] if k else [])
        samples[model] = [ids[i] for i in indices]
        report['models'][model] = {'source_path': str(source), 'numeric_entry_count': len(ids),
                                  'minimum_id': ids[0] if ids else None, 'maximum_id': ids[-1] if ids else None,
                                  'nonnumeric_entry_count': len(other), 'nonnumeric_names': sorted(other),
                                  'sample_ids': samples[model], 'census_complete': not any(e['source'] == str(source) for e in report['errors'])}
        (args.output / 'sample_ids.json').write_text(json.dumps(samples, indent=2) + '\n')
        (args.output / 'inventory.json').write_text(json.dumps(report, indent=2) + '\n')
        print(model, len(ids), 'numeric entries', flush=True)
    archive_dir = args.source / '50Am_data'
    try:
        with os.scandir(archive_dir) as entries:
            archive_paths = sorted(Path(entry.path) for entry in entries if entry.name.endswith('.tar.gz'))
        for path in archive_paths:
            stat = path.stat()
            report['archives'].append({'path': str(path), 'bytes': stat.st_size,
                                       'size_modulo_512': stat.st_size % 512,
                                       'mtime_utc': datetime.datetime.fromtimestamp(stat.st_mtime, datetime.timezone.utc).isoformat(),
                                       'integrity_verified': False})
    except OSError as exc:
        report['errors'].append({'source': str(archive_dir), 'error': str(exc)})
    report['archive_total_bytes'] = sum(a['bytes'] for a in report['archives'])
    log_findings = []
    excerpts = []
    pattern = re.compile(r'error|fail|cancel|time.?limit|truncat|unexpected|no space|complete|finished|exit|killed', re.I)
    for relative in ['logs', '50Am_data/logs', '50Am_data/slurm_log', '50Am_data/output']:
        directory = args.source / relative
        try:
            with os.scandir(directory) as entries:
                names = sorted(entry.name for entry in entries)
            for name in names[:96]:
                path = directory / name
                try:
                    with path.open('rb') as handle:
                        start = handle.read(4096)
                        size = handle.seek(0, os.SEEK_END)
                        if size > 4096:
                            handle.seek(max(4096, size - 4096))
                            end = handle.read(4096)
                        else:
                            end = b''
                    snippet = start.decode('utf-8', errors='replace')
                    if end:
                        snippet += '\n[TAIL; MIDDLE OMITTED IF FILE >8192 BYTES]\n' + end.decode('utf-8', errors='replace')
                    excerpts.append(f'FILE {path}\nBYTES {size}\n{snippet}\n')
                    log_findings.append({'path': str(path), 'bytes': size, 'fully_read': size <= 8192,
                                         'matching_lines': [line[:500] for line in snippet.splitlines() if pattern.search(line)][:30]})
                except OSError as exc:
                    log_findings.append({'path': str(path), 'error': str(exc)})
        except OSError as exc:
            report['errors'].append({'source': str(directory), 'error': str(exc)})
    report['log_sampling'] = {'scope': 'First 96 sorted direct entries per log folder; first and last 4096 bytes per file.', 'findings': log_findings}
    (args.output / 'archive_log_excerpts.txt').write_text('\n'.join(excerpts))
    (args.output / 'inventory.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'archive_count': len(report['archives']), 'archive_total_bytes': report['archive_total_bytes'], 'errors': report['errors']}), flush=True)


if __name__ == '__main__':
    main()
