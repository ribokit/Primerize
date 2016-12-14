from base import *

import os
import unittest


class TestPrimerize3D(unittest.TestCase):

    def test_default_normal(self):
        job_3d = prm_3d.design(INPUT['SEQ_P4P6'], structures=[INPUT['STR_P4P6']], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_single=False, is_fillWT=False, is_force=True)
        self.assertTrue(job_3d.is_success)
        self.assertDictEqual(job_3d._params, OUTPUT['3D']['normal']['param'])
        for i in range(job_3d.get('N_PRIMER')):
            self.assertEqual(len(job_3d._data['plates'][i][0]), OUTPUT['3D']['normal']['plate'][str(i + 1)])

    def test_default_diff(self):
        job_3d = prm_3d.design(INPUT['SEQ_P4P6'], structures=[INPUT['STR_P4P6_1'], INPUT['STR_P4P6_2']], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_single=True, is_fillWT=False, is_force=True)
        self.assertTrue(job_3d.is_success)
        for i in range(job_3d.get('N_PRIMER')):
            for j in range(job_3d.get('N_PLATE')):
                self.assertEqual(len(job_3d._data['plates'][i][j]), OUTPUT['3D']['diff']['plate'][str(j + 1)][str(i + 1)])

    def test_default_exclude(self):
        job_3d = prm_3d.design(INPUT['SEQ_P4P6'], structures=[INPUT['STR_P4P6_1'], INPUT['STR_P4P6_2']], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_exclude=True, is_single=True, is_fillWT=True, is_force=True)
        self.assertTrue(job_3d.is_success)
        for i in range(job_3d.get('N_PRIMER')):
            self.assertEqual(len(job_3d._data['plates'][i][0]), OUTPUT['3D']['exclude']['plate'])

    def test_str_mismatch(self):
        job_3d = prm_3d.design(INPUT['SEQ_P4P6'], structures=[INPUT['STR_P4P6_MISMATCH']], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_single=True, is_fillWT=True, is_force=True)
        self.assertTrue(job_3d.is_success)
        self.assertEqual(job_3d.get('warning'), [(89, 110), (90, 109), (93, 106), (94, 105), (97, 102), (113, 123), (114, 122), (85, 129), (86, 128), (75, 143), (76, 142), (77, 141), (80, 139), (81, 138), (82, 137), (64, 153), (67, 150), (68, 149), (69, 148), (60, 156), (55, 162), (58, 159), (175, 195), (176, 194), (177, 193), (178, 192), (179, 191), (181, 189), (168, 201), (169, 200), (171, 198), (163, 206), (164, 205)])
        job_3d.save('structures')
        os.remove('primer_structures.txt')


if __name__ == '__main__':
    unittest.main()
