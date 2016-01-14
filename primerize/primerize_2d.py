import argparse
import math
import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as pyplot
import time
import traceback

from util import *
from primerize_1d import Primerize_1D


class Design_2D(object):
    def __init__(self, sequence, name, is_success, primer_set, params, data):
        (self.sequence, self.name, self.is_success, self.primer_set, self._params, self._data) = (sequence, name, is_success, primer_set, params, data)

    def __repr__(self):
        pass

    def __str__(self):
        return self.echo()


    def get(self, keyword):
        keyword = keyword.upper()
        if self._params.has_key(keyword):
            return self._params[keyword]
        elif self._params.has_key(keyword.lower()):
            return self._params[keyword.lower()]
        elif keyword == 'PRIMER':
            return self._data['assembly'].primers
        elif keyword == 'CONSTRUCT':
            return self._data['construct_names']
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized keyword \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (keyword, self.__class__))


    def save(self, keyword='', path='./', name=None):
        if self.is_success:
            if name is None: name = self.name
            keyword = keyword.lower()
            if keyword == 'excel':
                save_plates_excel(self._data['plates'], self._params['N_PLATE'], self._params['N_PRIMER'], name, path)
            elif keyword == 'image':
                save_plate_layout(self._data['plates'], self._params['N_PLATE'], self._params['N_PRIMER'], name, path)
            elif keyword == 'plain' or keyword == 'text':
                save_construct_key(self._data['construct_names'], name, path)
            else:
                raise AttributeError('\033[41mERROR\033[0m: Unrecognized keyword \033[92m%s\033[0m for \033[94m%s.sav()\033[0m.\n' % (keyword, self.__class__)) 
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of keyword \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (keyword, self.__class__)) 


    def echo(self):
        if self.is_success:
            output = ''
            for i in xrange(len(self.plates[0])):
                for j in xrange(len(self.plates)):
                    output += 'Plate \033[95m%d\033[0m; Primer \033[92m%d\033[0m\n' % (i + 1, j + 1)
                    output += self.plates[j][i].print_constructs(self.primer_set[j])
            return output[:-1]
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of keyword \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (keyword, self.__class__)) 



class Primerize_2D(object):
    def __init__(self, offset=0, which_muts=[], which_libs=[1], prefix='lib'):
        self.prefix = prefix
        self.offset = offset
        self.which_muts = which_muts
        self.which_libs = which_libs

    def __repr__(self):
        return repr(self.__dict__)

    def __str__(self):
        return repr(self.__dict__)


    def get(self, keyword):
        keyword = keyword.lower()
        if hasattr(self, keyword):
            return getattr(self, keyword)
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized keyword \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (keyword, self.__class__))


    def set(self, keyword, value):
        keyword = keyword.lower()
        if hasattr(self, keyword):
            if keyword == 'prefix':
                self.prefix = str(value)
            elif keyword == 'offset' and isinstance(value, (int, float)):
                self.offset = int(value)
            elif keyword == 'which_libs' and (isinstance(value, list) and all(isinstance(x, (float, int)) for x in value) and all(x in (1, 2, 3) for x in value)):
                self.which_libs = sorted(set(value))
            elif keyword == 'which_muts' and (isinstance(value, list) and all(isinstance(x, (float, int)) for x in value)):
                self.which_muts = sorted(set(value))
            else:
                raise ValueError('\033[41mERROR\033[0m: Illegal value \033[95m%s\033[0m for keyword \033[92m%s\033[0m for \033[94m%s.set()\033[0m.\n' % (value, keyword, self.__class__)) 
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized keyword \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (keyword, self.__class__))


    def reset(self):
        self.prefix = 'lib'
        self.offset = 0
        self.which_muts = []
        self.which_libs = [1]


    def design(self, sequence, primer_set=None, offset=None, which_muts=None, which_libs=None, prefix=None):
        if primer_set is None: primer_set = self.primer_set
        if offset is None: offset = self.offset
        if which_muts is None: which_muts = self.which_muts
        if which_libs is None: which_libs = self.which_libs
        if prefix is None: prefix = self.prefix

        if len(primer_set) % 2:
            raise ValueError('\033[41mERROR\033[0m: Illegal length \033[95m%s\033[0m of value for params \033[92mprimer_set\033[0m for \033[94m%s.set()\033[0m.\n' % (len(primer_set), self.__class__)) 

        name = prefix
        sequence = RNA2DNA(sequence)
        N_BP = len(self.sequence)

        is_success = True
        assembly = {}
        for i in xrange(len(primer_set)):
            primer_set[i] = RNA2DNA(primer_set[i])
        if not primer_set:
            prm = Primerize_1D()
            res = prm.design(sequence)
            if res.is_success:
                primer_set = res.primer_set
            else:
                is_success = False
                print '\033[41mFAIL\033[0m: \033[31mNO Solution\033[0m found under given contraints.\n'
                params = {'offset': offset, 'which_muts': which_muts, 'which_libs': which_libs, 'N_BP': N_BP}
                data = {'plates': [], 'assembly': [], 'construct_names': []}
                return Design_2D(sequence, name, is_success, primer_set, params, data)
        if not which_muts:
            which_muts = range(1 - offset, N_BP + 1 - offset)
        construct_names = list(' ' * (len(which_muts) + 1))

        N_primers = len(primer_set)
        N_constructs = 1 + len(which_muts)
        N_plates = int(math.floor((N_constructs - 1) / 96.0) + 1)


        (primers, is_success) = get_primer_index(primer_set, sequence)
        if not is_success: 
            print '\033[41mFAIL\033[0m: \033[31mMismatch\033[0m of given \033[92mprimer_set\033[0m for given \033[92msequence\033[0m.\n'
            params = {'offset': offset, 'which_muts': which_muts, 'which_libs': which_libs, 'N_PRIMER': N_primers, 'N_PLATE': N_plates, 'N_CONSTRUCT': N_constructs, 'N_BP': N_BP}
            data = {'plates': [], 'assembly': [], 'construct_names': []}
            return Design_2D(sequence, name, is_success, primer_set, params, data)
        assembly = draw_assembly(sequence, primers, name)

        plates = [[Plate_96Well() for i in xrange(self.N_plates)] for i in xrange(self.N_primers)]
        print 'Filling out sequences ...'

        try:
            for p in xrange(N_primers):
                for l_pos in xrange(len(which_libs)):
                    # lib should be a number: 1, 2, or 3 for the different possible mutations.
                    lib = which_libs[l_pos]

                    for m_pos in xrange(-1, len(which_muts)):
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
                            plates[p][plate_num].set_well(well_tag, well_name, mut_primer)

            print '\033[92mSUCCESS\033[0m: Primerize 2D design() finished.\n'
        except:
            is_success = False
            print traceback.format_exc()
            print '\033[41mERROR\033[0m: Primerize 2D design() encountered error.\n'

        params = {'offset': offset, 'which_muts': which_muts, 'which_libs': which_libs, 'N_PRIMER': N_primers, 'N_PLATE': N_plates, 'N_CONSTRUCT': N_constructs, 'N_BP': N_BP}
        data = {'plates': plates, 'assembly': assembly, 'construct_names': construct_names}
        return Design_2D(sequence, name, is_success, primer_set, params, data)



class Plate_96Well(object):
    def __init__(self):
        self.coords = set()
        self.data = {}

    def set_well(self, coord, tag, primer):
        if coord_to_num(coord) == -1:
            print 'Invalid 96 well coordinate: %s.' % coord
        else:
            self.coords.add(coord)
            self.data[coord_to_num(coord)] = (tag, primer)

    def get_well(self, coord):
        if coord_to_num(coord) == -1:
            print 'Invalid 96 well coordinate: %s.' % coord
            return
        else:
            if coord in self.coords:
                return self.data[coord_to_num(coord)]
            else:
                return

    def get_count(self):
        return len(self.coords)

    def print_constructs(self, ref_primer=''):
        return print_primer_plate(self, ref_primer)

    def print_layout(self, file_name='./', title=''):
        fig = pyplot.figure()
        pyplot.axes().set_aspect('equal')
        pyplot.axis([0, 13.875, 0, 9.375])
        pyplot.xticks([x * 1.125 + 0.75 for x in range(12)], [str(x + 1) for x in range(12)], fontsize=14)
        pyplot.yticks([y * 1.125 + 0.75 for y in range(8)], list('ABCDEFGH'), fontsize=14)
        pyplot.suptitle(title, fontsize=16, fontweight='bold')
        ax = pyplot.gca()
        for edge in ('bottom', 'top', 'left', 'right'):
            ax.spines[edge].set_color('w')
        ax.invert_yaxis()
        ax.xaxis.set_ticks_position('top')
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False

        (x_green, x_violet, x_gray, y_green, y_violet, y_gray) = ([], [], [], [], [], [])
        for i in xrange(8):
            for j in xrange(12):
                num = i + j * 8 + 1
                if num_to_coord(num) in self.coords:
                    if 'WT' in self.data[num][0]:
                        x_green.append(j * 1.125 + 0.75)
                        y_green.append(i * 1.125 + 0.75)
                    else:
                        x_violet.append(j * 1.125 + 0.75)
                        y_violet.append(i * 1.125 + 0.75)
                else:
                    x_gray.append(j * 1.125 + 0.75)
                    y_gray.append(i * 1.125 + 0.75)
        pyplot.scatter(x_gray, y_gray, 961, c='#ffffff', edgecolor='#333333', linewidth=5)
        pyplot.scatter(x_violet, y_violet, 961, c='#ecddf4', edgecolor='#c28fdd', linewidth=5)
        pyplot.scatter(x_green, y_green, 961, c='#beebde', edgecolor='#29be92', linewidth=5)

        matplotlib.rcParams['svg.fonttype'] = 'none'
        pyplot.savefig(file_name, orientation='landscape', format='svg')



def design_primers_2D(sequence, primer_set=None, offset=None, which_muts=None, which_libs=None, prefix=None):
    prm = Primerize_2D()
    res = prm.design(sequence, primer_set, offset, which_muts, which_libs, prefix)
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
    group2.add_argument('-e', '--excel', action='store_true', dest='is_excel', help='Write Results to Excel File(s)')
    group2.add_argument('-i', '--image', action='store_true', dest='is_image', help='Save Layout to Image File(s)')
    group2.add_argument('-h', '--help', action='help', help='Show this Help Message')
    args = parser.parse_args()

    t0 = time.time()
    # if args.primer_set == None:
    #     args.primer_set = []
    # else:
    #     args.primer_set = args.primer_set[0]
    # if args.which_libs == None: args.which_libs = [1]
    # (which_muts, _, _) = get_mut_range(args.mut_start, args.mut_end, args.offset, args.sequence)

    res = design_primers_2D(args.sequence, args.primer_set, args.offset, which_muts, args.which_libs, args.prefix)
    if res.is_success:
        if not res.is_quiet:
            print res
        if args.is_excel:
            res.save('excel')
            res.save('text')
        if args.is_image:
            res.save('image')

    print 'Time elapsed: %.1f s.' % (time.time() - t0)


if __name__ == "__main__":
    main()

