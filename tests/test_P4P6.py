import primerize

import os
import simplejson
import unittest

prm_1d = primerize.Primerize_1D()
prm_2d = primerize.Primerize_2D()
prm_3d = primerize.Primerize_3D()

INPUT = simplejson.load(open(os.path.dirname(os.path.abspath(__file__)) + '/arg.json', 'r'))
OUTPUT = simplejson.load(open(os.path.dirname(os.path.abspath(__file__)) + '/res.json', 'r'))


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

        job_1d.save()
        os.remove('primer.txt')

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


class TestPrimerize2D(unittest.TestCase):

    def test_default(self):
        job_2d = prm_2d.design(INPUT['SEQ_P4P6'], primer_set=[], is_force=True)
        self.assertTrue(job_2d.is_success)
        self.assertDictEqual(job_2d._params, OUTPUT['2D']['default']['param'])
        for i in range(job_2d.get('N_PRIMER')):
            for j in range(job_2d.get('N_PLATE')):
                if str(i + 1) in OUTPUT['2D']['default']['plate'][str(j + 1)]:
                    self.assertEqual(len(job_2d._data['plates'][i][j]), OUTPUT['2D']['default']['plate'][str(j + 1)][str(i + 1)])
                else:
                    self.assertFalse(len(job_2d._data['plates'][i][j]))



# class TestPrimerize3D(unittest.TestCase):

#     def test_default_normal(self):
#         job_3d = prm_3d.design(INPUT['SEQ_P4P6'], primer_set=[], structures=[INPUT['STR_P4P6']], is_force=True)
#         self.assertTrue(job_3d.is_success)


if __name__ == '__main__':
    unittest.main()
