from base import *

import unittest


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

    def test_explicit(self):
        job_2d = prm_2d.design(INPUT['SEQ_P4P6'], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_force=True)
        self.assertTrue(job_2d.is_success)
        for i in range(job_2d.get('N_PRIMER')):
            for j in range(job_2d.get('N_PLATE')):
                if str(i + 1) in OUTPUT['2D']['explicit']['plate'][str(j + 1)]:
                    self.assertEqual(len(job_2d._data['plates'][i][j]), OUTPUT['2D']['explicit']['plate'][str(j + 1)][str(i + 1)])
                else:
                    self.assertFalse(len(job_2d._data['plates'][i][j]))


if __name__ == '__main__':
    unittest.main()
