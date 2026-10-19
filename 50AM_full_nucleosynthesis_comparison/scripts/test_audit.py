import unittest
import json,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
from common import canon,species,unique_rows
from population_audit import stats
import common,population_audit

class AuditTests(unittest.TestCase):
    def test_names(self):
        self.assertEqual(canon('X_{26Al}'),'26Al')
        self.assertEqual(species('X_{26Al}'),('isotope',13,26))
        self.assertEqual(species('X_{Al}'),('element',13,None))
        self.assertEqual(species('X_{n}'),('isotope',0,1))
    def test_ids_not_row_order_and_duplicate_exclusion(self):
        ids,rows=unique_rows(np.array([8,2,-1,9,2,4]))
        np.testing.assert_array_equal(ids,[4,8,9]);np.testing.assert_array_equal(rows,[5,0,3])
    def test_zero_missing_tiny_and_tolerance(self):
        old=np.array([0.,1.,np.nan,1e-310,2.,3.])
        new=np.array([1.,0.,2.,2e-310,2.+1e-6,3.])
        r,d=stats(old,new,np.arange(6),np.ones(6),np.ones(6),1e-5)
        self.assertEqual(r['newly_nonzero'],1);self.assertEqual(r['newly_zero'],1)
        self.assertEqual(r['missing_old_only'],1);self.assertEqual(r['exact'],1)
        self.assertEqual(r['approximate_not_exact'],1);self.assertEqual(r['significant'],3)
        self.assertTrue(np.isnan(d['ratio'][0]));self.assertTrue(np.isnan(d['log10_ratio'][1]))
        self.assertAlmostEqual(d['log10_ratio'][3],np.log10(2),places=11)
    def test_all_missing_yield_is_unavailable(self):
        r,d=stats(np.array([np.nan]),np.array([0.]),np.array([4]),np.ones(1),np.ones(1),1e-5)
        self.assertEqual(r['finite_pairs'],0);self.assertIsNone(r['new_yield_shared_msun'])
    def test_streaming_parser_preserves_malformed_and_missing(self):
        parent=common.ROOT/'tables/test_fixtures';parent.mkdir(exist_ok=True)
        root=Path(tempfile.mkdtemp(dir=parent));(root/'tables').mkdir();(root/'provenance').mkdir()
        source=root/'fixture.csv';source.write_text('id, X_{26Al}, X_{n}\n3, 1e-310, NA\ny2, 0.3y, 0\n1, 0, NA\n')
        with patch.object(common,'ROOT',root),patch.object(population_audit,'ROOT',root):population_audit.ingest('test',source)
        meta=json.loads((root/'provenance/test_scan.json').read_text());self.assertEqual(meta['rows'],3);self.assertEqual(meta['invalid_ids'],1)
        a=np.load(root/'tables/test_values.npy');self.assertEqual(a[0,1],1e-310);self.assertTrue(np.isnan(a[1,1]));self.assertEqual(a[1,2],0)
if __name__=='__main__':unittest.main()
