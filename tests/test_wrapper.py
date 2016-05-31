from base import *

import glob
import os
import unittest

class TestDesignSingle(unittest.TestCase):

    def setUp(self):
        self.job_1d = prm_1d.design(INPUT['SEQ_P4P6'], MIN_TM=INPUT['MIN_TM'], NUM_PRIMERS=INPUT['NUM_PRM'], MIN_LENGTH=INPUT['MIN_LEN'], MAX_LENGTH=INPUT['MAX_LEN'], prefix='primer')

    def test_save(self):
        self.job_1d.save()
        os.remove('primer.txt')

    def test_echo(self):
        self.job_1d.echo()

    def test_repr(self):
        print(repr(self.job_1d))

    def test_get(self):
        self.assertEqual(map(lambda x: list(x), self.job_1d.get('WARNING')), OUTPUT['1D']['warning'])
        self.assertListEqual(self.job_1d.get('MISPRIME'), OUTPUT['1D']['misprime'])
        self.assertEqual(self.job_1d.get('N_BP'), OUTPUT['1D']['param']['N_BP'])
        self.assertEqual(self.job_1d.get('MIN_TM'), OUTPUT['1D']['param']['MIN_TM'])
        self.assertEqual(self.job_1d.get('COL_SIZE'), OUTPUT['1D']['param']['COL_SIZE'])


class TestDesignPlate2(unittest.TestCase):

    def setUp(self):
        self.job_2d = prm_2d.design(INPUT['SEQ_P4P6'], primer_set=[], is_force=True)

    def test_save(self):
        self.job_2d.save()
        files = glob.glob("lib_*.xls")
        files.extend(glob.glob("lib_*.txt"))
        files.extend(glob.glob("Lib1_*.svg"))
        for f in files:
            os.remove(f)

    def test_echo(self):
        self.job_2d.echo()

    def test_repr(self):
        print(repr(self.job_2d))

    def test_get(self):
        self.assertEqual(self.job_2d.get('N_BP'), OUTPUT['2D']['default']['param']['N_BP'])
        self.assertEqual(self.job_2d.get('N_CONSTRUCT'), OUTPUT['2D']['default']['param']['N_CONSTRUCT'])
        self.assertEqual(self.job_2d.get('type'), OUTPUT['2D']['default']['param']['type'])


class TestDesignPlate3(unittest.TestCase):

    def setUp(self):
        self.job_3d = prm_3d.design(INPUT['SEQ_P4P6'], structures=[INPUT['STR_P4P6_1'], INPUT['STR_P4P6_2']], primer_set=INPUT['PRIMER_SET_P4P6'], offset=INPUT['OFFSET_P4P6'], which_muts=which_muts, which_lib=[int(INPUT['LIB_P4P6'])], prefix="primer", is_single=True, is_fillWT=True, is_force=True)

    def test_save(self):
        self.job_3d.save()
        files = glob.glob("primer_*.xls")
        files.extend(glob.glob("primer_*.txt"))
        files.extend(glob.glob("Lib1_*.svg"))
        for f in files:
            os.remove(f)

    def test_echo(self):
        self.job_3d.echo()

    def test_repr(self):
        print(repr(self.job_3d))

    def test_get(self):
        self.assertListEqual(self.job_3d.get('STRUCTURE'), [INPUT['STR_P4P6_1'], INPUT['STR_P4P6_2']])


if __name__ == '__main__':
    unittest.main()
