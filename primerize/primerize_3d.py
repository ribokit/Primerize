import argparse
import math
import time
import traceback

if __package__ is None or not __package__:
    from util import *
    from primerize_1d import Primerize_1D
    from wrapper import Design_Single, Design_Plate
else:
    from .util import *
    from .primerize_1d import Primerize_1D
    from .wrapper import Design_Single, Design_Plate


class Primerize_3D(object):
    def __init__(self, offset=0, N_mutations=1, which_lib=1, is_single=True, is_fillWT=False, COL_SIZE=142, prefix='lib'):
        self.prefix = prefix
        self.offset = offset
        self.N_mutations = max(min(N_mutations, 3), 1)
        self.which_lib = max(min(which_lib, 4), 1)
        self.is_single = bool(is_single)
        self.is_fillWT = bool(is_fillWT)
        self.COL_SIZE = max(COL_SIZE, 0)

    def __repr__(self):
        return repr(self.__dict__)

    def __str__(self):
        return repr(self.__dict__)


    def get(self, key):
        key = key.lower()
        if hasattr(self, key):
            return getattr(self, key)
        elif key == 'col_size':
            return self.COL_SIZE
        elif key == 'n_mutations':
            return self.N_mutations
        elif key == 'is_fillwt':
            return self.is_fillWT
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def set(self, key, value):
        key = key.lower()
        if hasattr(self, key):
            if key == 'prefix':
                self.prefix = str(value)
            elif key == 'offset' and isinstance(value, (int, float)):
                self.offset = int(value)
            elif key == 'n_mutations' and isinstance(value, (float, int)) and value in (1, 2, 3):
                self.N_mutations = int(value)
            elif key == 'which_lib' and isinstance(value, (float, int)) and value in (1, 4):
                self.which_lib = int(value)
            elif key == 'is_single':
                self.is_single = bool(value)
            elif key == 'is_fillwt':
                self.is_fillWT = bool(value)
            elif key == 'col_size' and isinstance(value, int) and value > 0:
                self.COL_SIZE = int(value)
            else:
                raise ValueError('\033[41mERROR\033[0m: Illegal value \033[95m%s\033[0m for key \033[92m%s\033[0m for \033[94m%s.set()\033[0m.\n' % (value, key, self.__class__))
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def reset(self):
        self.prefix = 'lib'
        self.offset = 0
        self.N_mutations = 1
        self.which_lib = 1
        self.is_single = True
        self.is_fillWT = False
        self.COL_SIZE = 142


    def design(self, sequence, primer_set=[], structures=[], offset=None, N_mutations=None, which_lib=None, which_muts=[], prefix=None, is_single=None, is_fillWT=False, is_force=False):
        if isinstance(sequence, Design_Single):
            design_1d = sequence
            sequence = design_1d.sequence
            primer_set = design_1d.primer_set
            prefix = design_1d.name

        offset = self.offset if offset is None else offset
        structures = [structures] if isinstance(structures, str) else structures
        N_mutations = self.N_mutations if N_mutations is None else N_mutations
        which_lib = self.which_lib if which_lib is None else which_lib
        is_single = self.is_single if is_single is None else is_single
        is_fillWT = self.is_fillWT if is_fillWT is None else is_fillWT
        prefix = self.prefix if prefix is None else prefix

        if len(primer_set) % 2:
            raise ValueError('\033[41mERROR\033[0m: Illegal length \033[95m%s\033[0m of value for params \033[92mprimer_set\033[0m for \033[94m%s.design()\033[0m.\n' % (len(primer_set), self.__class__))
        if not structures:
            raise ValueError('\033[41mERROR\033[0m: Missing input \033[92mstructures\033[0m for \033[94m%s.design()\033[0m.\n' % self.__class__)

        name = prefix
        sequence = RNA2DNA(sequence)
        N_BP = len(sequence)
        params = {'offset': offset, 'which_lib': which_lib, 'is_single': is_single, 'is_fillWT': is_fillWT, 'N_MUTATION': N_mutations, 'N_BP': N_BP, 'type': 'Mutation/Rescue'}
        data = {'plates': [], 'assembly': [], 'constructs': [], 'bps': []}

        is_success = True
        assembly = {}
        for i in range(len(primer_set)):
            primer_set[i] = RNA2DNA(primer_set[i])
        if not primer_set:
            if is_force:
                prm = Primerize_1D()
                res = prm.design(sequence)
                if res.is_success:
                    primer_set = res.primer_set
                else:
                    is_success = False
                    print('\033[41mFAIL\033[0m: \033[91mNO Solution\033[0m found under given contraints.\n')
            else:
                print('\033[93mWARNING\033[0m: Please run \033[34mPrimerize_1D.design()\033[0m first to get a solution for \033[92mprimer_set\033[0m.\n')
                is_success = False

        if not is_success:
            return Design_Plate({'sequence': sequence, 'name': name, 'is_success': is_success, 'primer_set': primer_set, 'structures': structures, 'params': params, 'data': data})

        if not which_muts:
            which_muts = list(range(1 - offset, N_BP + 1 - offset))
        else:
            which_muts = [x for x in which_muts if x >= 1 - offset and x < N_BP + 1 - offset]
        which_lib = which_lib[0] if isinstance(which_lib, list) else which_lib
        N_primers = len(primer_set)
        params.update({'which_muts': which_muts, 'which_lib': which_lib, 'N_PRIMER': N_primers})

        (primers, is_success) = get_primer_index(primer_set, sequence)
        if not is_success:
            print('\033[41mFAIL\033[0m: \033[91mMismatch\033[0m of given \033[92mprimer_set\033[0m for given \033[92msequence\033[0m.\n')
            return Design_Plate({'sequence': sequence, 'name': name, 'is_success': is_success, 'primer_set': primer_set, 'structures': structures, 'params': params, 'data': data})

        assembly = Assembly(sequence, primers, name, self.COL_SIZE)
        constructs = Construct_List()
        data.update({'assembly': assembly, 'constructs': constructs})

        bps = diff_bps(structures)
        for pair in list(bps):
            if not (pair[0] - offset in which_muts and pair[1] - offset in which_muts):
                bps.remove(pair)
        if not bps:
            print('\033[41mFAIL\033[0m: \033[91mNo\033[0m base-pairs exist within given \033[92mstructures\033[0m and \033[92mwhich_muts\033[0m.\n')
            return Design_Plate({'sequence': sequence, 'name': name, 'is_success': False, 'primer_set': primer_set, 'structures': structures, 'params': params, 'data': data})


        N_constructs = (len(bps) - N_mutations + 1) * (is_single * 2 + 1) + 1
        constructs.push('WT')
        N_plates = int(math.floor((N_constructs - 1) / 96.0) + 1)
        plates = [[Plate_96Well(which_lib) for i in range(N_plates)] for j in range(N_primers)]

        for i in range(len(bps) - N_mutations + 1):
            (mut_list_l, mut_list_r) = ([], [])

            for j in range(N_mutations):
                if sequence[bps[i + j][0] - 1] == 'G' and sequence[bps[i + j][1] - 1] == 'T':
                    mut_list_l.append('G%dC' % (bps[i + j][0] - offset))
                    mut_list_r.append('T%dG' % (bps[i + j][1] - offset))
                elif sequence[bps[i + j][0] - 1] == 'T' and sequence[bps[i + j][1] - 1] == 'G':
                    mut_list_l.append('T%dG' % (bps[i + j][0] - offset))
                    mut_list_r.append('G%dC' % (bps[i + j][1] - offset))
                else:
                    mut_list_l.append('%s%d%s' % (sequence[bps[i + j][0] - 1], bps[i + j][0] - offset, get_mutation(sequence[bps[i + j][0] - 1], which_lib)))
                    mut_list_r.append('%s%d%s' % (sequence[bps[i + j][1] - 1], bps[i + j][1] - offset, get_mutation(sequence[bps[i + j][1] - 1], which_lib)))

            if is_single:
                constructs.push(mut_list_l)
                constructs.push(mut_list_r)
            constructs.push(mut_list_l + mut_list_r)

        try:
            plates = mutate_primers(plates, primers, primer_set, offset, constructs, which_lib, is_fillWT)
            print('\033[92mSUCCESS\033[0m: Primerize 3D design() finished.\n')
        except:
            is_success = False
            print(traceback.format_exc())
            print('\033[41mERROR\033[0m: Primerize 3D design() encountered error.\n')

        params.update({'N_PLATE': N_plates, 'N_CONSTRUCT': N_constructs})
        data.update({'plates': plates, 'constructs': constructs, 'bps': bps})
        return Design_Plate({'sequence': sequence, 'name': name, 'is_success': is_success, 'primer_set': primer_set, 'structures': structures, 'params': params, 'data': data})



def design_primers_3D(sequence, primer_set=[], offset=None, structures=[], N_mutations=None, which_muts=None, which_lib=None, is_single=None, is_fillWT=None, prefix=None):
    prm = Primerize_3D()
    res = prm.design(sequence, primer_set, offset, structures, N_mutations, which_lib, which_muts, prefix, is_single, is_fillWT, True)
    return res


def main():
    parser = argparse.ArgumentParser(description='\033[92mPrimerize 3D Mutation/Rescue Plate Design\033[0m', epilog='\033[94mby Siqi Tian, 2016\033[0m', add_help=False)
    parser.add_argument('sequence', type=str, help='DNA Template Sequence')
    parser.add_argument('-p', metavar='prefix', type=str, help='Display Name of Construct', dest='prefix', default='lib')
    group1 = parser.add_argument_group('advanced options')
    group1.add_argument('-s', metavar='PRIMER_SET', type=str, nargs='+', help='Set of Primers for Assembly (Default runs Primerize 1D)', dest='primer_set', action='append')
    group1.add_argument('-r', metavar='STRUCTURES', type=str, nargs='+', help='Set of Structures for Mutation/Rescue, difference between each pair of structures are included, use "quote"', dest='structures', action='append')
    group1.add_argument('-o', metavar='OFFSET', type=int, help='Sequence Numbering Offset', dest='offset', default=0)
    group1.add_argument('-l', metavar='MUT_START', type=int, help='First Position of Mutagenesis (Inclusive), numbering with OFFSET applied', dest='mut_start', default=None)
    group1.add_argument('-u', metavar='MUT_END', type=int, help='Last Position of Mutagenesis (Inclusive), numbering with OFFSET applied', dest='mut_end', default=None)
    group1.add_argument('-n', metavar='N_MUTATION', type=int, choices=(1, 2, 3), help='Number of mutations for "single mutant" {1, 2, 3}', dest='N_mutations', default=1)
    group1.add_argument('-w', metavar='LIB', type=int, choices=(1, 4), help='Mutation Library Choices {1, 4}, GU pairs will be mutated to CG pairs', dest='which_lib', default=1)
    group1.add_argument('-g', '--single', action='store_true', dest='is_single', help='Include single mutants')
    group1.add_argument('-f', '--fill', action='store_true', dest='is_fillWT', help='Fill Wild-Type primers')
    group2 = parser.add_argument_group('commandline options')
    group2.add_argument('-q', '--quiet', action='store_true', dest='is_quiet', help='Suppress Results Printing to stdout')
    group2.add_argument('-e', '--excel', action='store_true', dest='is_excel', help='Write Order Table to Excel File(s)')
    group2.add_argument('-i', '--image', action='store_true', dest='is_image', help='Save Layout to Image File(s)')
    group2.add_argument('-t', '--text', action='store_true', dest='is_text', help='Save Construct and Assembly to Text File(s)')
    group2.add_argument('-h', '--help', action='help', help='Show this Help Message')
    args = parser.parse_args()

    t0 = time.time()
    args.primer_set = [] if args.primer_set is None else args.primer_set[0]
    (which_muts, _, _) = get_mut_range(args.mut_start, args.mut_end, args.offset, args.sequence)
    args.structures = [] if args.structures is None else args.structures[0]


    res = design_primers_3D(args.sequence, args.primer_set, args.structures, args.offset, args.N_mutations, which_muts, args.which_lib, args.is_single, args.is_fillWT, args.prefix)
    if res.is_success:
        if not args.is_quiet:
            print(res)
        if args.is_excel:
            res.save('table')
        if args.is_image:
            res.save('image')
        if args.is_text:
            res.save('construct')
            res.save('assembly')

    print('Time elapsed: %.1f s.' % (time.time() - t0))


if __name__ == "__main__":
    main()

