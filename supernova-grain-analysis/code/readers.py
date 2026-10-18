"""Strict readers for the observed SNSPH and PRISM formats; sources are read-only."""
from pathlib import Path
import re
import numpy as np

SYMBOLS = ('n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn '
           'Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd '
           'In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu '
           'Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu '
           'Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()
Z_BY_SYMBOL = {s: i for i, s in enumerate(SYMBOLS)}
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
TIME_HEADER = re.compile(r'timestep\s+(\d+)\s+time\(s\)\s*('+NUMBER+r')\s+'
                         r'temperature\(GK\)\s*('+NUMBER+r')\s+density\(g/cm\^3\)\s*('+NUMBER+r')\s*$')

def fnum(s):
    return float(s.replace('D', 'E').replace('d', 'e'))

def sdf_header(path):
    """Read declarations; reject ambiguous schema, byte order or record count."""
    path = Path(path)
    with path.open('rb') as f:
        buf = f.read(65536)
    marker = re.search(br'# SDF-EOH[ \t]*\r?\n', buf)
    if not buf.startswith(b'# SDF') or marker is None:
        raise ValueError(f'No bounded SDF header: {path}')
    offset = marker.end()
    header = buf[:offset].decode('ascii')
    clean = re.sub(r'/\*.*?\*/', '', header, flags=re.S)
    match = re.search(r'struct\s*\{(.*?)\}\s*\[(\d+)\]\s*;', clean, re.S)
    if not match:
        raise ValueError('Unsupported SDF struct declaration')
    bo = re.search(r'byteorder\s*=\s*(0x[\da-fA-F]+)', clean)
    endian = {'0x78563412': '<', '0x12345678': '>'}.get(bo.group(1) if bo else '')
    if endian is None:
        raise ValueError('Unknown byte order')
    types = {'double':'f8', 'float':'f4', 'int':'i4', 'unsigned int':'u4'}
    fields = []
    for declaration in match.group(1).split(';'):
        if not declaration.strip():
            continue
        d = re.fullmatch(r'\s*(unsigned\s+int|double|float|int)\s+(.+?)\s*', declaration, re.S)
        if not d:
            raise ValueError(f'Unsupported declaration: {declaration}')
        typ = types[' '.join(d.group(1).split())]
        for name in d.group(2).split(','):
            name = name.strip()
            if not re.fullmatch(r'\w+', name):
                raise ValueError(f'Unsupported field: {name}')
            fields.append((name, endian+typ))
    dtype = np.dtype(fields, align=False)
    meta = {k: fnum(v) for k,v in re.findall(r'(?:int|float|double)\s+(\w+)\s*=\s*('+NUMBER+r')\s*;', clean[:match.start()])}
    n = int(match.group(2))
    if 'npart' in meta and int(meta['npart']) != n:
        raise ValueError('npart disagrees with struct count')
    actual = path.stat().st_size
    expected = offset + n * dtype.itemsize
    if actual != expected:
        raise ValueError(f'Payload length mismatch: actual={actual}, expected={expected}')
    return {'path':str(path), 'header':header, 'offset':offset, 'npart':n,
            'record_bytes':dtype.itemsize, 'dtype':dtype, 'meta':meta, 'file_bytes':actual}

def open_sdf(path):
    info = sdf_header(path)
    return info, np.memmap(path, dtype=info['dtype'], mode='r', offset=info['offset'], shape=(info['npart'],))

def units(meta):
    """Return factors plus provenance. Defaults are evidence-backed, still inferred."""
    if 'lenCF' in meta and 'lengthCF' in meta and meta['lenCF'] != meta['lengthCF']:
        raise ValueError('Conflicting length factors')
    length = meta.get('lengthCF', meta.get('lenCF', 6.955e10))
    mass = meta.get('massCF', 1.9889e27)
    time = meta.get('timeCF', 100.)
    if min(length, mass, time) <= 0:
        raise ValueError('Nonpositive unit conversion')
    explicit = 'massCF' in meta and 'timeCF' in meta and ('lenCF' in meta or 'lengthCF' in meta)
    return {'length_cm':length, 'mass_g':mass, 'time_s':time,
            'density_g_cm3':mass/length**3, 'velocity_cm_s':length/time,
            'acceleration_cm_s2':length/time**2,
            'provenance':'header' if explicit else 'inferred_from_shared_code_defaults'}

def table(text, keycols):
    result = {}
    for line_no, line in enumerate(text.splitlines(), 1):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        p = line.split()
        if len(p) != keycols + 1:
            raise ValueError(f'Malformed abundance row {line_no}: {line[:100]}')
        key = tuple(int(s) for s in p[:keycols])
        if key in result:
            raise ValueError(f'Duplicate abundance key {key}')
        value = fnum(p[-1])
        if not np.isfinite(value) or value < 0:
            raise ValueError(f'Invalid abundance at {key}: {value}')
        if keycols == 2 and (key[0] < 0 or key[1] < max(1,key[0])):
            raise ValueError(f'Invalid nuclear Z,A: {key}')
        result[key] = value
    return result

def read_extent(path):
    out = set()
    seen = set()
    for row in Path(path).read_text().splitlines():
        if not row.strip() or row.lstrip().startswith('#'):
            continue
        z, lo, hi = map(int, row.split())
        if z in seen or z < 0 or lo < 0 or hi < lo:
            raise ValueError(f'Invalid/duplicate extent row: {row}')
        seen.add(z)
        out.update((z,z+n) for n in range(lo,hi+1))
    return out

def time_blocks(lines):
    """Yield header and isotope dictionary for each block without keeping history."""
    current, values = None, {}
    for raw in lines:
        line = raw.decode('ascii') if isinstance(raw, bytes) else raw
        m = TIME_HEADER.search(line)
        if m:
            if current is not None:
                yield current, values
            current = {'step':int(m[1]), 'time_s':fnum(m[2]),
                       'temp_GK':fnum(m[3]), 'rho_g_cm3_as_written':fnum(m[4])}
            if not all(np.isfinite(v) for v in current.values()) or current['temp_GK']<=0 or current['rho_g_cm3_as_written']<=0:
                raise ValueError(f'Invalid thermodynamic header: {line[:150]}')
            values = {}
        elif line.lstrip().startswith('timestep'):
            raise ValueError(f'Unrecognized time header: {line[:150]}')
        elif current is not None and line.strip():
            row = table(line, 2)
            if values.keys() & row.keys():
                raise ValueError('Duplicate species within time block')
            values.update(row)
    if current is not None:
        yield current, values

def request_kind(label):
    if label == 'X_{n}':
        return 'isotope', (0,1)
    m = re.fullmatch(r'(\d+)([A-Z][a-z]?)', label)
    if m:
        return 'isotope', (Z_BY_SYMBOL[m[2]],int(m[1]))
    if label in Z_BY_SYMBOL:
        return 'element', (Z_BY_SYMBOL[label],)
    return 'physical', None
