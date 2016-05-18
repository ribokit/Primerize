from base import *

import unittest


class TestThermo(unittest.TestCase):

    def setUp(self):
        self.NN = primerize.thermo.Nearest_Neighbor_Parameter()

    def test_singleton(self):
        self.assertEqual(id(self.NN), id(primerize.Nearest_Neighbor))

    def test_NN(self):
        self.assertEqual(self.NN.T, 273.15 + 37)
        self.assertEqual(self.NN.delS_init, -5.7)
        self.assertListEqual(self.NN.delH_AT_closing_penalty.tolist(), [2.2, 0.0, 0.0, 2.2])

    def test_calc_Tm(self):
        self.assertEqual(round(primerize.thermo.calc_Tm('ACCAATTATCATCAAGTATT'), 2), 65.49)
        self.assertEqual(round(primerize.thermo.calc_Tm('TATATTATTTTTATTTATAT'), 2), 51.80)
        self.assertEqual(round(primerize.thermo.calc_Tm('CGAGCGGCAGAGAGCGCGGT', DNA_conc=0.2e-6, monovalent_conc=0.1, divalent_conc=0.0015), 2), 76.91)


if __name__ == '__main__':
    unittest.main()
