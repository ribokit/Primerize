from base import *

import unittest


class TestPrimerize1D(unittest.TestCase):

    def test_default(self):
        job_1d = prm_1d.design(INPUT['SEQ_P4P6'])
        self.assertTrue(job_1d.is_success)
        self.assertListEqual(job_1d.primer_set, OUTPUT['1D']['primer'])
        self.assertListEqual(map(lambda x: round(x, 2), job_1d._data['assembly'].Tm_overlaps), OUTPUT['1D']['Tm'])
        for i, coord in enumerate(OUTPUT['1D']['coord']):
            self.assertListEqual(job_1d._data['assembly'].primers[i, :].tolist(), coord)
        self.assertListEqual(map(lambda x: list(x), job_1d._data['warnings']), OUTPUT['1D']['warning'])
        self.assertDictEqual(job_1d._params, OUTPUT['1D']['param'])

    def test_default_explicit(self):
        job_1d = prm_1d.design(INPUT['SEQ_P4P6'], MIN_TM=INPUT['MIN_TM'], NUM_PRIMERS=INPUT['NUM_PRM'], MIN_LENGTH=INPUT['MIN_LEN'], MAX_LENGTH=INPUT['MAX_LEN'], prefix='primer')
        self.assertTrue(job_1d.is_success)
        self.assertListEqual(job_1d.primer_set, OUTPUT['1D']['primer'])

    def test_Tm_65(self):
        job_1d = prm_1d.design(INPUT['SEQ_P4P6'], MIN_TM=65)
        self.assertTrue(job_1d.is_success)
        self.assertListEqual(job_1d.primer_set, INPUT['PRIMER_SET_P4P6'])

    def test_Tm_70(self):
        job_1d = prm_1d.design(INPUT['SEQ_P4P6'], MIN_TM=70)
        self.assertFalse(job_1d.is_success)


if __name__ == '__main__':
    unittest.main()
