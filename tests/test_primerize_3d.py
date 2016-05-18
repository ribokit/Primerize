from base import *

import unittest


class TestPrimerize3D(unittest.TestCase):

    def test_default_normal(self):
        job_3d = prm_3d.design(INPUT['SEQ_P4P6'], structures=[INPUT['STR_P4P6']], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_single=False, is_fillWT=False, is_force=True)
        self.assertTrue(job_3d.is_success)
        self.assertDictEqual(job_3d._params, OUTPUT['3D']['normal']['param'])
        for i in range(job_3d.get('N_PRIMER')):
            self.assertEqual(len(job_3d._data['plates'][i][0]), OUTPUT['3D']['normal']['plate'][str(i + 1)])

    def test_default_diff(self):
        job_3d = prm_3d.design(INPUT['SEQ_P4P6'], structures=[INPUT['STR_P4P6_1'], INPUT['STR_P4P6_2']], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_single=True, is_fillWT=True, is_force=True)
        self.assertTrue(job_3d.is_success)
        for i in range(job_3d.get('N_PRIMER')):
            self.assertEqual(len(job_3d._data['plates'][i][0]), OUTPUT['3D']['diff']['plate'])


if __name__ == '__main__':
    unittest.main()
