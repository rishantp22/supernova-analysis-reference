"""Package the delivered table and notes for transfer without changing sources."""
import hashlib
import json
from pathlib import Path
import time
import zipfile

ROOT = Path('/scratch/rpraka15')
WORK = Path(__file__).resolve().parent
META = WORK / 'plotting_format'
DEST = ROOT / '50Am_UCF_OneDrive.zip'
FILES = [(ROOT / '50Am_plotting_complete.out', '50Am_AM_export/50Am_plotting_complete.out')]
for name in ['complete_model_README.txt', 'README.txt', 'column_map.csv',
             'complete_model_validation.json', 'complete_model_coverage.csv',
             'complete_model_preview.out', 'validation.json']:
    FILES.append((META / name, '50Am_AM_export/notes/' + name))
FILES.append((WORK / 'mass_fraction_findings.txt', '50Am_AM_export/notes/mass_fraction_findings.txt'))


def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True)


def main():
    if DEST.exists(): raise FileExistsError(DEST)
    temporary = DEST.with_suffix('.zip.partial')
    manifest = []
    with zipfile.ZipFile(temporary, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=3) as archive:
        for path, name in FILES:
            before = path.stat()
            digest = hashlib.sha256()
            log(f'Packaging {path.name}')
            with path.open('rb') as source, archive.open(name, 'w', force_zip64=True) as target:
                for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
                    target.write(block); digest.update(block)
            after = path.stat()
            assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
            manifest.append({'archive_path':name, 'source_path':str(path),
                             'bytes':before.st_size, 'sha256':digest.hexdigest()})
        expected = json.loads((META / 'complete_model_validation.json').read_text())['output_sha256']
        assert manifest[0]['sha256'] == expected
        archive.writestr('50Am_AM_export/manifest.json', json.dumps(manifest, indent=2))
    log('Verifying all archived files against their SHA-256 checksums')
    with zipfile.ZipFile(temporary) as archive:
        for record in manifest:
            digest = hashlib.sha256()
            with archive.open(record['archive_path']) as source:
                for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
                    digest.update(block)
            assert digest.hexdigest() == record['sha256']
    temporary.rename(DEST)
    print(json.dumps({'archive':str(DEST), 'bytes':DEST.stat().st_size,
                      'files_including_manifest':len(manifest)+1,
                      'all_file_sha256_checks_passed':True,
                      'saved_to_OneDrive':False}, indent=2), flush=True)


if __name__ == '__main__':
    main()
