from base import *

import unittest


class TestThermo(unittest.TestCase):

    def test_singleton(self):
        NN = primerize.thermo.Nearest_Neighbor_Wrapper()
        self.assertEqual(id(NN), id(primerize.Nearest_Neighbor))


if __name__ == '__main__':
    unittest.main()
