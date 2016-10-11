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


if __name__ == '__main__':
    unittest.main()
