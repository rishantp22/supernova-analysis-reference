"""Small fixtures for correctness of the streamed, name/ID-based comparison."""
import tempfile
import unittest
from pathlib import Path

import compare_outputs as audit


class ComparisonTests(unittest.TestCase):
    def test_ids_columns_duplicates_units_and_isotope_difference(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            h = ['Id','x','y','z','vx','vy','vz','ax','ay','az','mass','h',
                 'Density','Peak temp','Peak density','Y_{e}','X_{n}','H','1H','4He']
            values = ['1']*15 + ['NA','NA','0.1','0.2']
            (p/'a').write_text('\t'.join(h)+'\n'+''.join(
                str(i)+'\t'+'\t'.join(values)+'\n' for i in [2,1,2,4]))
            names = [audit.canonical(s) for s in h]
            bv = ['NA' if s == 'NA' else str(float(s)*audit.FACTORS.get(n,1))
                  for n,s in zip(names[1:],values)]
            bv[-1] = '0.3'
            hb = ['id']+[('X_{'+s+'}') if s in ['H','1H','4He'] else s for s in h[1:]]
            order = [0]+list(range(len(h)-1,0,-1))
            rows = []
            for i in [1,3,2]:
                v = [str(i)]+bv
                rows.append(', '.join(v[j] for j in order))
            (p/'b').write_text(', '.join(hb[j] for j in order)+'\n'+'\n'.join(rows)+'\n')
            a,ap = audit.scan(p/'a',b'\t',p,'a')
            b,bp = audit.scan(p/'b',b',',p,'b')
            result,columns = audit.compare(a,b,ap,bp,p)
            self.assertEqual((a['duplicate_ids'],a['duplicate_extra_rows']),(1,1))
            self.assertEqual((result['only_final_chunks_ids'],result['only_plotting_ids']),(1,1))
            self.assertEqual((result['compared_unique_ids'],result['ambiguous_shared_duplicate_ids']),(1,1))
            self.assertEqual(result['isotope_exact_numeric_mismatches'],1)
            self.assertEqual(result['all_column_missing_mismatches'],0)
            self.assertFalse(result['column_order_equal_after_name_mapping'])
            self.assertEqual(next(c for c in columns if c['column']=='x')['converted_mismatches_rtol_1e9'],0)

    def test_rejects_invalid_numbers_and_nonfinite_values(self):
        for token in [b'bad', b'inf', b'nan']:
            with self.subTest(token=token), self.assertRaises((ValueError,Warning)):
                audit.numbers([b'1\t'+token],2)

    def test_missing_is_not_zero_and_tiny_isotopes_are_retained(self):
        values = audit.numbers([b'NA\t0\t1e-315'],3)
        self.assertTrue(audit.np.isnan(values[0,0]))
        self.assertEqual(values[0,1],0)
        self.assertNotEqual(values[0,2],0)

    def test_invalid_tokens_are_separate_from_missing(self):
        values, invalid = audit.audited_numbers([b'NA\t-0.045y2421\t0'],3)
        self.assertEqual(invalid.tolist(),[[False,True,False]])
        self.assertTrue(audit.np.isnan(values[0,0]))
        self.assertTrue(audit.np.isnan(values[0,1]))
        self.assertEqual(values[0,2],0)


if __name__ == '__main__':
    unittest.main()
