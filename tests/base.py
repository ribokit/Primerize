import primerize

import os
import simplejson

__all__ = ['primerize', 'prm_1d', 'prm_2d', 'prm_3d', 'prm_cm', 'INPUT', 'OUTPUT', 'which_muts']

prm_1d = primerize.Primerize_1D
prm_2d = primerize.Primerize_2D
prm_3d = primerize.Primerize_3D
prm_cm = primerize.Primerize_Custom

INPUT = simplejson.load(open(os.path.dirname(os.path.abspath(__file__)) + '/input.json', 'r'))
OUTPUT = simplejson.load(open(os.path.dirname(os.path.abspath(__file__)) + '/output.json', 'r'))
(which_muts, _, _) = primerize.util.get_mut_range(INPUT['MIN_MUTS_P4P6'], INPUT['MAX_MUTS_P4P6'], INPUT['OFFSET_P4P6'], INPUT['SEQ_P4P6'])
