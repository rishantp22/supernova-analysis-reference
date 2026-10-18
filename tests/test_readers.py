"""Small format-contract tests; never access master simulation data."""
import io
from pathlib import Path
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from readers import open_sdf, sdf_header, table, time_blocks, read_extent, request_kind


class ReaderTests(unittest.TestCase):
    def setUp(self):
        # Preserve even test fixtures: the project owner requested no deletions.
        self.root = Path(tempfile.mkdtemp(prefix='supernova-reader-tests-'))

    def make_sdf(self, endian, count=1, declared_count=None):
        marker = '0x78563412' if endian == '<' else '0x12345678'
        header = ('# SDF 1.0\nparameter byteorder = ' + marker + ';\n'
                  'int npart = ' + str(count if declared_count is None else declared_count) + ';\n'
                  'double tpos = 2.5;\n'
                  'struct { double x; float mass; double temp; unsigned int ident; } [' + str(count) + '];\n'
                  '# SDF-EOH  \n').encode('ascii')
        # Deliberately packed: the second double is at byte 12, not byte 16.
        payload = struct.pack(endian + 'dfdI', -3.25, 0.125, 9.75e8, 4000000001)
        path = self.root / ('little.sdf' if endian == '<' else 'big.sdf')
        path.write_bytes(header + payload * count)
        return path

    def test_both_endians_packed_unsigned_and_readonly(self):
        for endian in ['<', '>']:
            with self.subTest(endian=endian):
                path = self.make_sdf(endian)
                info, records = open_sdf(path)
                self.assertEqual(info['record_bytes'], 24)
                self.assertEqual(info['dtype'].fields['temp'][1], 12)
                self.assertEqual(info['npart'], 1)
                self.assertEqual(info['meta']['tpos'], 2.5)
                self.assertEqual(float(records['x'][0]), -3.25)
                self.assertEqual(float(records['temp'][0]), 9.75e8)
                self.assertEqual(int(records['ident'][0]), 4000000001)
                self.assertFalse(records.flags.writeable)
                del records

    def test_truncated_payload_and_disagreeing_count_rejected(self):
        path = self.make_sdf('<')
        path.write_bytes(path.read_bytes()[:-1])
        with self.assertRaisesRegex(ValueError, 'Payload length mismatch'):
            sdf_header(path)
        path = self.make_sdf('<', declared_count=2)
        with self.assertRaisesRegex(ValueError, 'npart disagrees'):
            sdf_header(path)

    def test_duplicate_or_nonfinite_abundance_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            table('6 12 0.1\n6 12 0.2\n', 2)
        for value in ['nan', 'inf', '-0.1']:
            with self.subTest(value=value), self.assertRaises(ValueError):
                table('6 12 ' + value, 2)

    def test_missing_sparse_species_is_not_zero_or_outside_extent(self):
        extent_path = self.root / 'extent.dat'
        extent_path.write_text('0 1 1\n6 2 19\n')
        extent = read_extent(extent_path)
        values = table('6 12 1.0D-2\n', 2)
        self.assertEqual(request_kind('13C'), ('isotope', (6, 13)))
        self.assertEqual(request_kind('X_{n}'), ('isotope', (0, 1)))
        self.assertEqual(request_kind('C'), ('element', (6,)))
        self.assertEqual(request_kind('Y_{e}'), ('physical', None))
        self.assertEqual(values[(6, 12)], 0.01)
        self.assertNotIn((6, 13), values)
        self.assertIsNone(values.get((6, 13)))
        self.assertIn((6, 13), extent)
        self.assertNotIn((6, 26), extent)
        self.assertIn((6, 8), extent)
        self.assertIn((6, 25), extent)

    def test_glued_time_labels_and_duplicate_block_species(self):
        header = ('timestep      2   time(s) 0.100000D-14   '
                  'temperature(GK)0.672383E+01   density(g/cm^3)0.238415E+07\n')
        blocks = list(time_blocks(io.StringIO('  Z A X(Z,A)\n' + header + '6 12 1e-2\n')))
        self.assertEqual(len(blocks), 1)
        head, values = blocks[0]
        self.assertEqual(head['step'], 2)
        self.assertAlmostEqual(head['temp_GK'], 6.72383)
        self.assertAlmostEqual(head['rho_g_cm3_as_written'], 2384150)
        self.assertEqual(values[(6, 12)], 0.01)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            list(time_blocks(io.StringIO(header + '6 12 0.1\n6 12 0.2\n')))

    def test_invalid_time_headers_and_extent_rejected(self):
        for density in ['1e999', '-1', '1e+', '0']:
            with self.subTest(density=density), self.assertRaises(ValueError):
                list(time_blocks(io.StringIO('timestep 2 time(s) 1 temperature(GK) 2 density(g/cm^3)'+density+'\n')))
        path=self.root/'extent.dat'
        for text in ['6 19 2\n','6 -1 3\n','6 2 19\n6 2 20\n']:
            path.write_text(text)
            with self.assertRaises(ValueError):read_extent(path)


if __name__ == '__main__':
    unittest.main()
