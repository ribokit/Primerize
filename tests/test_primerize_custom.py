from base import *

import unittest


class TestPrimerizeCustom(unittest.TestCase):

    def test_default_normal(self):
        mut_list = primerize.util.Construct_List()
        mut_list.push(['T120C'])
        mut_list.push(['G119A'])
        mut_list.push(['G119A', 'T120C'])
        mut_list.push('WT')
        job_cm = prm_cm.design(INPUT['SEQ_P4P6'], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], mut_list=mut_list, prefix="primer", is_force=True)
        self.assertTrue(job_cm.is_success)
        print(repr(job_cm))
        print(job_cm.echo())
        self.assertDictEqual(job_cm._params, OUTPUT['custom']['normal']['param'])

    def test_default_mismatch(self):
        mut_list = primerize.util.Construct_List()
        mut_list.push(['A120C'])
        job_cm = prm_cm.design(INPUT['SEQ_P4P6'], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], mut_list=mut_list, prefix="primer", is_force=True)
        self.assertFalse(job_cm.is_success)

    def test_default_outbound(self):
        mut_list = primerize.util.Construct_List()
        mut_list.push(['A12C'])
        mut_list.push(['A11C'])
        mut_list.push(['A10C'])
        mut_list.push(['C9A'])
        job_cm = prm_cm.design(INPUT['SEQ_P4P6'], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], mut_list=mut_list, prefix="primer", is_force=True)
        self.assertTrue(job_cm.is_success)
        self.assertEqual(job_cm.get('N_CONSTRUCT'), 5)
        self.assertEqual(len(job_cm._data['plates'][0][0]), 1)

    def test_merge(self):
        job_2d = prm_2d.design(INPUT['SEQ_P4P6'], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_force=True)
        job_3d = prm_3d.design(INPUT['SEQ_P4P6'], structures=[INPUT['STR_P4P6']], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_single=True, is_fillWT=False, is_force=True)
        mut_list = job_2d.get('CONSTRUCT')
        repeated = mut_list.merge(job_3d.get('CONSTRUCT'))
        self.assertEqual(len(repeated), 87)

        job_cm = prm_cm.design(INPUT['SEQ_P4P6'], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], mut_list=mut_list)
        self.assertTrue(job_cm.is_success)
        print(repr(job_cm))
        print(job_cm.echo())
        self.assertDictEqual(job_cm._params, OUTPUT['custom']['merge']['param'])
        for i in range(job_cm.get('N_PRIMER')):
            for j in range(job_cm.get('N_PLATE')):
                if str(i + 1) in OUTPUT['custom']['merge']['plate'][str(j + 1)]:
                    self.assertEqual(len(job_cm._data['plates'][i][j]), OUTPUT['custom']['merge']['plate'][str(j + 1)][str(i + 1)])
                else:
                    self.assertFalse(len(job_cm._data['plates'][i][j]))


if __name__ == '__main__':
    unittest.main()
