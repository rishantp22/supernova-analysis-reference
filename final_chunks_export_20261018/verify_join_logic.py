"""Small adversarial checks for ID/species joins, independent of production data."""
import argparse
import importlib.util
import json
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('exporter', HERE / 'build_final_chunks.py')
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)
columns = (exporter.ROOT / 'List of isotopes_elements_variables.txt').read_text().splitlines()
species = [key for name in columns for kind, key in [exporter.request_kind(name)]
           if kind == 'isotope' and name != 'X_{n}']
rows = []
for ident, mass in [(123, 4.), (456, 7.)]:
    for category_index, category in enumerate(exporter.CATEGORIES):
        carbon = (0.2 if ident == 123 else 0.4) + category_index * 0.1
        for z, a in species:
            value = carbon if (z, a) == (6, 12) else (1-carbon if (z, a) == (8, 16) else 0.)
            rows.append([ident, category, category_index + 1, float(category_index + 1), z, a, value])
frame_a = pd.DataFrame(rows, columns=['ID','category','actual_timestep','actual_time','A','B','C_sum'])
frame_a.category = frame_a.category.astype(exporter.CAT_TYPE)
frame_b = frame_a[['ID','category','A','B']].copy()
frame_b['mass_times_Csum'] = frame_a.C_sum * frame_a.ID.map({123:4.,456:7.})
physical = pd.DataFrame({name: ['1','2'] for name in exporter.PHYSICAL.values()}, index=[123,456])
physical['ID'] = ['123','456']
physical['Mass'] = ['4','7']
peaks = pd.DataFrame({'Peak temp':[100.,200.], 'Peak density':[3.,4.],
                      'hydro_snapshot_count':[61,61]}, index=[123,456])
args = argparse.Namespace(category='last', elements='na')
results = []

for test in ['shuffled_keys', 'duplicate_key', 'missing_key']:
    a = frame_a.sample(frac=1, random_state=7).reset_index(drop=True)
    b = frame_b.sample(frac=1, random_state=13).reset_index(drop=True)
    if test == 'duplicate_key': b = pd.concat([b, b.iloc[[0]]], ignore_index=True)
    if test == 'missing_key': b = b.iloc[1:].copy()
    with tempfile.TemporaryDirectory(prefix='am_join_check_') as temp:
        exporter.WORK = Path(temp)
        (exporter.WORK / 'parts').mkdir()
        exporter.load_csv = lambda path, **kwargs: ((b if path.name.endswith('_B.csv') else a).copy(), {})
        try:
            record, ids = exporter.process_chunk(1, args, columns, physical, peaks)
        except ValueError as error:
            expected = 'Duplicate' if test == 'duplicate_key' else 'key mismatch'
            assert test != 'shuffled_keys' and expected in str(error), str(error)
        else:
            assert test == 'shuffled_keys'
            result = pd.read_csv(exporter.WORK / 'parts/0001.out', sep='\t').set_index('Id')
            assert list(result.index) == [123,456]
            np.testing.assert_allclose(result['12C'], [0.4,0.6])
            np.testing.assert_allclose(result['16O'], [0.6,0.4])
            np.testing.assert_allclose(result['mass'], [4.,7.])
            assert result['C'].isna().all() and result['X_{n}'].isna().all()
            assert record['selected_B_mass_mismatches'] == 0
            assert not record['A_B_same_row_order']
        results.append({'test':test, 'result':'PASS'})

with (HERE / 'join_logic_checks.json').open('x') as stream:
    json.dump(results, stream, indent=2)
print(json.dumps(results))
