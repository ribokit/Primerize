import argparse
import math
import time
import traceback

if __package__ is None or not __package__:
    from util import *
    from primerize_1d import Primerize_1D, Design_1D
else:
    from .util import *
    from .primerize_1d import Primerize_1D, Design_1D


class Design_2D(object):
    def __init__(self, sequence, name, is_success, primer_set, params, data):
        (self.sequence, self.name, self.is_success, self.primer_set, self._params, self._data) = (sequence, name, is_success, primer_set, params, data)

    def __repr__(self):
        return '\033[94m%s\033[0m {\n\033[95msequence\033[0m = \'%s\', \n\033[95mname\033[0m = \'%s\', \n\033[95mis_success\033[0m = \033[41m%s\033[0m, \n\033[95mprimer_set\033[0m = %s, \n\033[95mparams\033[0m = %s, \n\033[95mdata\033[0m = {\n    \033[92m\'construct_names\'\033[0m: \033[31mlist\033[0m(\033[31mstring\033[0m * %d), \n    \033[92m\'assembly\'\033[0m: %s, \n    \033[92m\'plates\'\033[0m: %s\n}' % (self.__class__, self.sequence, self.name, self.is_success, repr(self.primer_set), repr(self._params), len(self._data['construct_names']), repr(self._data['assembly']), repr(self._data['plates']))

    def __str__(self):
        return self.echo()


    def get(self, key):
        key = key.upper()
        if key in self._params:
            return self._params[key]
        elif key.lower() in self._params:
            return self._params[key.lower()]
        elif key == 'PRIMER':
            return self._data['assembly'].primers
        elif key == 'CONSTRUCT':
            return self._data['construct_names']
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def save(self, key='', path='./', name=None):
        if self.is_success:
            if name is None: name = self.name
            key = key.lower()
            if key == 'table':
                save_plates_excel(self._data['plates'], self._params['N_PLATE'], self._params['N_PRIMER'], name, path)
            elif key == 'image':
                save_plate_layout(self._data['plates'], self._params['N_PLATE'], self._params['N_PRIMER'], name, path)
            elif key == 'construct':
                save_construct_key(self._data['construct_names'], name, path)
            elif key == 'assembly':
                self._data['assembly'].save(path, name)

            elif not key:
                for key in ['table', 'image', 'construct', 'assembly']:
                    self.save(key, path, name)
            else:
                raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.save()\033[0m.\n' % (key, self.__class__))
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of key \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (key, self.__class__))


    def echo(self, key=''):
        if self.is_success:
            key = key.lower()
            if key == 'plate':
                output = ''
                for i in range(len(self._data['plates'][0])):
                    for j in range(len(self._data['plates'])):
                        output += 'Plate \033[95m%d\033[0m; Primer \033[92m%d\033[0m\n' % (i + 1, j + 1)
                        output += self._data['plates'][j][i].echo(self.primer_set[j])
                return output[:-1]
            elif key == 'assembly':
                return self._data['assembly'].echo()

            elif not key:
                return self.echo('assembly') + '\n\n' + self.echo('plate')
            else:
                raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.echo()\033[0m.\n' % (key, self.__class__))
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of key \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (key, self.__class__))



class Primerize_2D(object):
    def __init__(self, offset=0, which_muts=[], which_libs=[1], COL_SIZE=142, prefix='lib'):
        self.prefix = prefix
        self.offset = offset
        self.which_muts = which_muts
        self.which_libs = which_libs
        self.COL_SIZE = COL_SIZE

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
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def set(self, key, value):
        key = key.lower()
        if hasattr(self, key):
            if key == 'prefix':
                self.prefix = str(value)
            elif key == 'offset' and isinstance(value, (int, float)):
                self.offset = int(value)
            elif key == 'which_libs' and (isinstance(value, list) and all(isinstance(x, (float, int)) for x in value) and all(x in (1, 2, 3) for x in value)):
                self.which_libs = sorted(set(value))
            elif key == 'which_muts' and (isinstance(value, list) and all(isinstance(x, (float, int)) for x in value)):
                self.which_muts = sorted(set(value))
            elif key == 'col_size' and isinstance(value, int) and value > 0:
                self.COL_SIZE = int(value)
            else:
                raise ValueError('\033[41mERROR\033[0m: Illegal value \033[95m%s\033[0m for key \033[92m%s\033[0m for \033[94m%s.set()\033[0m.\n' % (value, key, self.__class__))
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def reset(self):
        self.prefix = 'lib'
        self.offset = 0
        self.which_muts = []
        self.which_libs = [1]
        self.COL_SIZE = 142


    def design(self, sequence, primer_set=[], offset=None, which_muts=None, which_libs=None, prefix=None, is_force=False):

        if isinstance(sequence, Design_1D):
            design_1d = sequence
            sequence = design_1d.sequence
            primer_set = design_1d.primer_set
            prefix = design_1d.name

        if offset is None: offset = self.offset
        if which_muts is None: which_muts = self.which_muts
        if which_libs is None: which_libs = self.which_libs
        if prefix is None: prefix = self.prefix

        if len(primer_set) % 2:
            raise ValueError('\033[41mERROR\033[0m: Illegal length \033[95m%s\033[0m of value for params \033[92mprimer_set\033[0m for \033[94m%s.set()\033[0m.\n' % (len(primer_set), self.__class__))

        name = prefix
        sequence = RNA2DNA(sequence)
        N_BP = len(sequence)

        is_success = True
        assembly = {}
        for i in range(len(primer_set)):
            primer_set[i] = RNA2DNA(primer_set[i])
        if not primer_set and is_force:
            prm = Primerize_1D()
            res = prm.design(sequence)
            if res.is_success:
                primer_set = res.primer_set
            else:
                is_success = False
                print('\033[41mFAIL\033[0m: \033[31mNO Solution\033[0m found under given contraints.\n')
        else:
            print('\033[93mWARNING\033[0m: Please run \033[31mPrimerize_1D.design()\033[0m first to get a solution for \033[94mprimer_set\033[0m.\n')
            is_success = False

        if not is_success:
            params = {'offset': offset, 'which_muts': which_muts, 'which_libs': which_libs, 'N_BP': N_BP}
            data = {'plates': [], 'assembly': [], 'construct_names': []}
            return Design_2D(sequence, name, is_success, primer_set, params, data)

        if not which_muts:
            which_muts = list(range(1 - offset, N_BP + 1 - offset))
        construct_names = list(' ' * (len(which_muts) + 1))

        N_primers = len(primer_set)
        N_constructs = 1 + len(which_muts)
        N_plates = int(math.floor((N_constructs - 1) / 96.0) + 1)


        (primers, is_success) = get_primer_index(primer_set, sequence)
        if not is_success:
            print('\033[41mFAIL\033[0m: \033[31mMismatch\033[0m of given \033[92mprimer_set\033[0m for given \033[92msequence\033[0m.\n')
            params = {'offset': offset, 'which_muts': which_muts, 'which_libs': which_libs, 'N_PRIMER': N_primers, 'N_PLATE': N_plates, 'N_CONSTRUCT': N_constructs, 'N_BP': N_BP}
            data = {'plates': [], 'assembly': [], 'construct_names': []}
            return Design_2D(sequence, name, is_success, primer_set, params, data)

        assembly = Assembly(sequence, primers, name, self.COL_SIZE)
        plates = [[Plate_96Well() for i in range(N_plates)] for i in range(N_primers)]
        print('Filling out sequences ...')

        try:
            for p in range(N_primers):
                for l_pos in range(len(which_libs)):
                    # lib should be a number: 1, 2, or 3 for the different possible mutations.
                    lib = which_libs[l_pos]

                    for m_pos in range(-1, len(which_muts)):
                        # which construct is this?
                        n = m_pos + 1
                        plate_num = int(math.floor(n / 96.0))
                        plate_pos = n % 96 + 1
                        well_tag = num_to_coord(plate_pos)

                        # m is actual position along sequence
                        if m_pos == -1:
                            m = -1
                        else:
                            m = offset + which_muts[m_pos] - 1

                        if (m >= primers[0, p] and m <= primers[1, p]) or m == -1:
                            wt_primer = primer_set[p]
                            mut_primer = wt_primer
                            if m == -1:
                                well_name = 'Lib%d-WT' % lib
                            else:
                                if primers[2, p] == -1:
                                    wt_primer = reverse_complement(wt_primer)
                                    mut_primer = reverse_complement(mut_primer)

                                m_shift = int(m - primers[0, p])
                                mut_primer = list(mut_primer)
                                mut_primer[m_shift] = get_mutation(wt_primer[m_shift], lib)
                                mut_primer = ''.join(mut_primer)

                                # Name, e.g., "C75A".
                                well_name = 'Lib%d-%s%d%s' % (lib, wt_primer[m_shift], which_muts[m_pos], mut_primer[m_shift])

                                if primers[2, p] == -1:
                                    wt_primer = reverse_complement(wt_primer)
                                    mut_primer = reverse_complement(mut_primer)

                            construct_names[n] = well_name
                            plates[p][plate_num].set(well_tag, well_name, mut_primer)

            print('\033[92mSUCCESS\033[0m: Primerize 2D design() finished.\n')
        except:
            is_success = False
            print(traceback.format_exc())
            print('\033[41mERROR\033[0m: Primerize 2D design() encountered error.\n')

        params = {'offset': offset, 'which_muts': which_muts, 'which_libs': which_libs, 'N_PRIMER': N_primers, 'N_PLATE': N_plates, 'N_CONSTRUCT': N_constructs, 'N_BP': N_BP}
        data = {'plates': plates, 'assembly': assembly, 'construct_names': construct_names}
        return Design_2D(sequence, name, is_success, primer_set, params, data)



def design_primers_2D(sequence, primer_set=[], offset=None, which_muts=None, which_libs=None, prefix=None):
    prm = Primerize_2D()
    res = prm.design(sequence, primer_set, offset, which_muts, which_libs, prefix, True)
    return res


def main():
    parser = argparse.ArgumentParser(description='\033[92mPrimerize 2D Mutate-and-Map Plate Design\033[0m', epilog='\033[94mby Siqi Tian, 2016\033[0m', add_help=False)
    parser.add_argument('sequence', type=str, help='DNA Template Sequence')
    parser.add_argument('-p', metavar='prefix', type=str, help='Display Name of Construct', dest='prefix', default='lib')
    group1 = parser.add_argument_group('advanced options')
    group1.add_argument('-s', metavar='PRIMER_SET', type=str, nargs='+', help='Set of Primers for Assembly (Default runs Primerize 1D)', dest='primer_set', action='append')
    group1.add_argument('-o', metavar='OFFSET', type=int, help='Sequence Numbering Offset', dest='offset', default=0)
    group1.add_argument('-l', metavar='MUT_START', type=int, help='First Position of Mutagenesis (Inclusive)', dest='mut_start', default=None)
    group1.add_argument('-u', metavar='MUT_END', type=int, help='Last Position of Mutagenesis (Inclusive)', dest='mut_end', default=None)
    group1.add_argument('-w', metavar='LIB', type=int, choices=(1, 2, 3), nargs='+', help='Mutation Library Choices {1, 2, 3}', dest='which_libs', action='append')
    group2 = parser.add_argument_group('commandline options')
    group2.add_argument('-q', '--quiet', action='store_true', dest='is_quiet', help='Suppress Results Printing to stdout')
    group2.add_argument('-e', '--excel', action='store_true', dest='is_excel', help='Write Order Table to Excel File(s)')
    group2.add_argument('-i', '--image', action='store_true', dest='is_image', help='Save Layout to Image File(s)')
    group2.add_argument('-t', '--text', action='store_true', dest='is_text', help='Save Construct and Assembly to Text File(s)')
    group2.add_argument('-h', '--help', action='help', help='Show this Help Message')
    args = parser.parse_args()

    t0 = time.time()
    if args.primer_set == None:
        args.primer_set = []
    else:
        args.primer_set = args.primer_set[0]
    (which_muts, _, _) = get_mut_range(args.mut_start, args.mut_end, args.offset, args.sequence)


    res = design_primers_2D(args.sequence, args.primer_set, args.offset, which_muts, args.which_libs, args.prefix)
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

