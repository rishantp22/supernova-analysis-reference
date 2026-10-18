"""Independent final sum/spot checks; originals and published export are read-only."""
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path('/scratch/rpraka15')
WORK = Path(__file__).resolve().parent
summary = json.loads((WORK / 'validation.json').read_text())
epochs = pd.read_csv(WORK / 'grain_epochs.csv')
assert len(epochs) == summary['rows'] and epochs.Id.is_unique
assert epochs.hydro_matched.all()
assert (epochs.hydro_snapshot_count == summary['peak_snapshot_count']).all()
minimum = epochs.loc[epochs.selected_isotope_sum.idxmin()]
maximum = epochs.loc[epochs.selected_isotope_sum.idxmax()]
targets = {int(minimum.Id): int(minimum.source_chunk), int(maximum.Id): int(maximum.source_chunk)}
for idx in [0, len(epochs)//2, len(epochs)-1]:
    row = epochs.iloc[idx]
    targets[int(row.Id)] = int(row.source_chunk)
symbols = ('n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn '
           'Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd '
           'In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu '
           'Hf Ta W Re Os Ir Pt').split()
source_values = {}
for ident, chunk in targets.items():
    values = {}
    with (ROOT / 'final_chunks' / f'chunk_{chunk}.csv').open() as handle:
        for row in csv.DictReader(handle):
            if int(row['ID']) == ident and row['category'] == 'last':
                isotope = row['B'] + symbols[int(row['A'])]
                assert isotope not in values
                values[isotope] = float(row['C_sum'])
    assert len(values) == 222
    source_values[ident] = values

checks = []
with (ROOT / 'final_chunks.out').open() as handle:
    fields = handle.readline().rstrip('\n').split('\t')
    expected_header = (ROOT / 'List of isotopes_elements_variables.txt').read_text().splitlines()
    assert fields == expected_header
    for line in handle:
        ident = int(line.partition('\t')[0])
        if ident not in source_values:
            continue
        row = dict(zip(fields, line.rstrip('\n').split('\t')))
        values = source_values[ident]
        for isotope, expected in values.items():
            assert float(row[isotope]) == expected, (ident, isotope, row[isotope], expected)
        checks.append({'Id': ident, 'source_chunk': targets[ident],
                       'isotope_fields_exactly_matched': len(values),
                       'independent_sum': math.fsum(values.values())})
assert len(checks) == len(targets)

digest = hashlib.sha256()
with (ROOT / 'final_chunks.out').open('rb') as handle:
    for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
        digest.update(block)
assert digest.hexdigest() == summary['output_sha256']

extreme = epochs[(epochs.selected_isotope_sum > 1.0001) | (epochs.selected_isotope_sum < 0.9999)]
extreme.to_csv(WORK / 'last_fraction_sums_outside_1e-4.csv', index=False, float_format='%.17g', mode='x')
quantiles = epochs.selected_isotope_sum.quantile([0., .01, .05, .5, .95, .99, 1.])
result = {'output_sha256_verified': True, 'rows_in_epoch_file': len(epochs),
          'all_export_IDs_have_physical_fields': True,
          'all_export_IDs_have_all_61_hydro_snapshots': True,
          'independent_source_spot_checks': checks,
          'last_sum_quantiles': {str(k):float(v) for k,v in quantiles.items()},
          'last_sums_above_1_plus_1e_4': epochs.loc[epochs.selected_isotope_sum > 1.0001,
              ['Id','source_chunk','actual_timestep','actual_time_as_written','selected_isotope_sum']].to_dict('records')}
with (WORK / 'final_verification.json').open('x') as handle:
    json.dump(result, handle, indent=2)

last = summary['category_totals']['last']
lines = [
    'Credibility checks for final_chunks.out', '',
    f'Output: {summary["rows"]:,} distinct IDs, {summary["columns"]} columns, last category.',
    'All 938 pairs join exactly by ID, category and isotope. No duplicates.',
    'Every exported ID has AM physical fields and all 61 saved hydro snapshots.',
    'All 624,194,514 B rows across the three categories match AM mass * C_sum within 2e-5 relative tolerance.',
    'The positive B/code-mass/X ratios are within approximately 5e-6 of one.', '',
    'UNITS',
    'Raw B values are in native mass units. The local SNSPH convention is 1 code mass = 1e-6 solar masses.',
    'Evidence: sn_utils.py:58 defines this convention; the original part2 script multiplies mass * C_sum;',
    'the original part3 aggregation script:251 applies the factor 1e-6 afterwards.',
    'AM SDF headers omit massCF, so this conversion is supported by local source code, not an explicit header.', '',
    'SUMS (DO NOT ADD THE THREE EPOCHS TOGETHER)',
]
for category, stats in summary['category_totals'].items():
    lines.append(f'{category}: B sum converted = {stats["sum_B_solar_masses"]:.12g} solar masses; '
                 f'per-ID selected-isotope sums min/max = {stats["sum_X_min"]:.12g} / {stats["sum_X_max"]:.12g}.')
lines += ['',
    f'For last, {last["sum_X_eq_one_within_1e_4"]:,}/{last["IDs"]:,} IDs are within 1e-4 of one.',
    f'{last["sum_X_below_1_minus_1e_4"]:,} are below 0.9999 and {last["sum_X_above_1_plus_1e_4"]} is above 1.0001.',
    f'The high outlier is ID {int(maximum.Id)} in chunk {int(maximum.source_chunk)}: sum = {maximum.selected_isotope_sum:.16g}.',
    'This was independently re-summed directly from the original chunk and preserved unchanged in the export.',
    'The selected isotope sums are not uniformly one. No values were renormalized.',
    f'The same {last["IDs"]:,} particles have total mass {last["sum_AM_particle_mass_solar"]:.12g} solar masses.',
    f'The difference from their last B sum is {last["sum_AM_particle_mass_solar"] - last["sum_B_solar_masses"]:.12g} solar masses.',
    f'The hydro snapshot contains {summary["physical_IDs_without_chunks"]:,} additional IDs absent from the chunks;',
    f'the full snapshot particle mass is {summary["physical_verification"]["all_hydro_mass_code"] * 1e-6:.12g} solar masses.',
    'These are population/species accounting sums, not independent proof of the total unbound/ejected mass.', '',
    'Elemental totals and free-neutron fraction are NA because the provided chunks cannot recover them.',
    'Physical fields retain native units. Peak fields are sampled maxima over the 61 AM hydro snapshots.',
    'Original source zeros remain zeros, including zeros inserted by the historical generator for missing species.',
    'All original files and the published export were left unchanged during these checks.',
    'See README.txt, mass_fraction_validation.json, normalization_by_id.csv and final_verification.json.'
]
with (WORK / 'mass_fraction_findings.txt').open('x') as handle:
    handle.write('\n'.join(lines) + '\n')
print(json.dumps(result, indent=2))
