"""Create and verify a new tar.gz containing both delivered AM tables."""
import hashlib
import io
import json
import os
from pathlib import Path
import tarfile
import time

from package_for_onedrive import FILES, META, ROOT, log

DEST = ROOT / '50Am_UCF_OneDrive.tar.gz'
FILES = list(FILES)
FILES.insert(1, (ROOT / '50Am_plotting_all_isotopes.out',
                 '50Am_AM_export/50Am_plotting_all_isotopes.out'))
BLOCK = 4 * 1024 * 1024


class HashingReader:
    def __init__(self, source):
        self.source = source
        self.digest = hashlib.sha256()

    def read(self, size=-1):
        data = self.source.read(size)
        self.digest.update(data)
        return data


def metadata(path):
    stat = path.stat()
    return (stat.st_size, stat.st_mtime_ns, stat.st_ino, stat.st_dev)


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    partial = Path(str(DEST) + '.partial')
    initial = {path: metadata(path) for path, _ in FILES}
    expected = {}
    for name in ('complete_model_validation.json', 'validation.json'):
        report = json.loads((META / name).read_text())
        expected[Path(report['output'])] = report['output_sha256']
    manifest = []
    with partial.open('xb') as target:
        with tarfile.open(fileobj=target, mode='w:gz', compresslevel=3,
                          format=tarfile.PAX_FORMAT, copybufsize=BLOCK) as archive:
            for path, name in FILES:
                log(f'Packaging {path.name}')
                info = archive.gettarinfo(str(path), arcname=name)
                if not info.isfile():
                    raise ValueError(f'Expected a regular file: {path}')
                with path.open('rb') as source:
                    reader = HashingReader(source)
                    archive.addfile(info, reader)
                checksum = reader.digest.hexdigest()
                if metadata(path) != initial[path]:
                    raise RuntimeError(f'Source changed during packaging: {path}')
                if path in expected and checksum != expected[path]:
                    raise RuntimeError(f'Previously validated checksum mismatch: {path}')
                manifest.append({'archive_path': name, 'source_path': str(path),
                                 'bytes': initial[path][0], 'sha256': checksum})
            payload = json.dumps(manifest, indent=2).encode() + b'\n'
            info = tarfile.TarInfo('50Am_AM_export/manifest.json')
            info.size = len(payload)
            info.mtime = int(time.time())
            info.mode = 0o644
            archive.addfile(info, io.BytesIO(payload))
        target.flush()
        os.fsync(target.fileno())

    log('Reading the complete archive back and checking every file checksum')
    records = {record['archive_path']: record for record in manifest}
    seen = set()
    manifest_verified = False
    with tarfile.open(partial, mode='r|gz') as archive:
        for member in archive:
            if not member.isfile() or member.name in seen:
                raise RuntimeError(f'Unexpected archive entry: {member.name}')
            seen.add(member.name)
            with archive.extractfile(member) as source:
                if member.name == '50Am_AM_export/manifest.json':
                    if source.read() != payload:
                        raise RuntimeError('Archive manifest mismatch')
                    manifest_verified = True
                    continue
                record = records[member.name]
                if member.size != record['bytes']:
                    raise RuntimeError(f'Archive size mismatch: {member.name}')
                digest = hashlib.sha256()
                for block in iter(lambda: source.read(BLOCK), b''):
                    digest.update(block)
                if digest.hexdigest() != record['sha256']:
                    raise RuntimeError(f'Archive checksum mismatch: {member.name}')
                log(f'Verified {member.name}')
    if not manifest_verified or seen != set(records) | {'50Am_AM_export/manifest.json'}:
        raise RuntimeError('Archive members did not match the manifest')
    if any(metadata(path) != initial[path] for path, _ in FILES):
        raise RuntimeError('Source metadata changed')
    # Publish without overwriting any existing archive.
    os.link(partial, DEST)
    partial.unlink()
    print(json.dumps({'archive': str(DEST), 'bytes': DEST.stat().st_size,
                      'files_including_manifest': len(seen),
                      'all_file_sha256_checks_passed': True,
                      'source_files_unchanged': True,
                      'saved_to_OneDrive': False}, indent=2), flush=True)


if __name__ == '__main__':
    main()
