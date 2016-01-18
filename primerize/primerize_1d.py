import argparse
import math
from numba import jit
import numpy
import os
import time
import traceback

if __package__ is None or not __package__:
    from misprime import *
    from thermo import *
    from util import *
else:
    from .misprime import *
    from .thermo import *
    from .util import *


class Design_1D(object):
    def __init__(self, sequence, name, is_success, primer_set, params, data):
        (self.sequence, self.name, self.is_success, self.primer_set, self._params, self._data) = (sequence, name, is_success, primer_set, params, data)

    def __repr__(self):
        return '\033[94m%s\033[0m {\n\033[95msequence\033[0m = \'%s\', \n\033[95mname\033[0m = \'%s\', \n\033[95mis_success\033[0m = \033[41m%s\033[0m, \n\033[95mprimer_set\033[0m = %s, \n\033[95mparams\033[0m = %s, \n\033[95mdata\033[0m = {\n    \033[92m\'misprime_score\'\033[0m: %s, \n    \033[92m\'assembly\'\033[0m: %s, \n    \033[92m\'warnings\'\033[0m: %s\n}' % (self.__class__, self.sequence, self.name, self.is_success, repr(self.primer_set), repr(self._params), repr(self._data['misprime_score']), repr(self._data['assembly']), repr(self._data['warnings']))

    def __str__(self):
        return self.echo()


    def get(self, key):
        key = key.upper()
        if key in self._params:
            return self._params[key]
        elif key == 'WARNING':
            return self._data['warnings']
        elif key == 'PRIMER':
            return self._data['asssembly'].primers
        elif key == 'MISPRIME':
            return self._data['misprime_score']
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__)) 


    def save(self, path='./', name=None):
        if self.is_success:
            if name is None: name = self.name
            f = open(os.path.join(path, '%s.txt' % name), 'w')

            f.write('Primerize Result\n\nINPUT\n=====\n%s\n' % self.sequence)
            f.write('#\nMIN_TM: %.1f\n' % self._params['MIN_TM'])
            if not self._params['NUM_PRIMERS']:
                f.write('NUM_PRIMERS: auto (unspecified)')
            else:
                f.write('NUM_PRIMERS: %d' % self._params['NUM_PRIMERS'])
            f.write('\nMAX_LENGTH: %d\nMIN_LENGTH: %d\n' % (self._params['MAX_LENGTH'], self._params['MIN_LENGTH']))

            f.write('\n\nOUTPUT\n======\n')
            lines = str(self).replace('\033[0m', '').replace('\033[100m', '').replace('\033[92m', '').replace('\033[93m', '').replace('\033[94m', '').replace('\033[95m', '').replace('\033[96m', '').replace('\033[41m', '')
            f.write(lines)
            f.write('#\n\n------/* IDT USER: for primer ordering, copy and paste to Bulk Input */------\n------/* START */------\n')
            for i in range(len(self.primer_set)):
                suffix = 'FR'[i % 2]
                f.write('%s-%d%s\t%s\t\t25nm\tSTD\n' % (self.name, i + 1, suffix, self.primer_set[i]))
            f.write('------/* END */------\n------/* NOTE: use "Lab Ready" for "Normalization" */------\n')
            f.close()
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of key \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (key, self.__class__)) 


    def echo(self, key=''):
        if self.is_success:
            key = key.lower()
            if key == 'misprime':
                output = ''
                for i in range(int(math.floor(self._params['N_BP'] / self._params['COL_SIZE'])) + 1):
                    output += '%s\n\033[92m%s\033[0m\n%s\n\n' % (self._data['misprime_score'][0][i * self._params['COL_SIZE']:(i + 1) * self._params['COL_SIZE']], self.sequence[i * self._params['COL_SIZE']:(i + 1) * self._params['COL_SIZE']], self._data['misprime_score'][1][i * self._params['COL_SIZE']:(i + 1) * self._params['COL_SIZE']])
                return output[:-1]
 
            elif key == 'warning':
                output = ''
                for i in range(len(self._data['warnings'])):
                    warning = self._data['warnings'][i]
                    p_1 = '\033[100m%d\033[0m%s' % (warning[0], primer_suffix(warning[0] - 1))
                    p_2 = ', '.join('\033[100m%d\033[0m%s' % (x, primer_suffix(x - 1)) for x in warning[3])
                    output += '\033[93mWARNING\033[0m: Primer %s can misprime with %d-residue overlap to position %s, which is covered by primers: %s\n' % (p_1.rjust(4), warning[1], str(int(warning[2])).rjust(3), p_2)
                return output[:-1]
 
            elif key == 'primer':
                output = '%s%s\tSEQUENCE\n' % ('PRIMERS'.ljust(20), 'LENGTH'.ljust(10))
                for i in range(len(self.primer_set)):
                    name = '%s-\033[100m%s\033[0m%s' % (self.name, i + 1, primer_suffix(i))
                    output += '%s%s\t%s\n' % (name.ljust(39), str(len(self.primer_set[i])).ljust(10), self.primer_set[i])
                return output[:-1]

            elif key == 'assembly':
                return self._data['assembly'].echo()
            elif not key:
                return self.echo('misprime') + '\n' + self.echo('assembly') + '\n' + self.echo('primer') + '\n\n' + self.echo('warning') + '\n'

            else:
                raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.echo()\033[0m.\n' % (key, self.__class__)) 
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of key \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (key, self.__class__)) 



class Primerize_1D(object):
    def __init__(self, MIN_TM=60.0, NUM_PRIMERS=0, MIN_LENGTH=15, MAX_LENGTH=60, COL_SIZE=142, WARN_CUTOFF=3, prefix='primer'):
        self.prefix = prefix
        self.MIN_TM = max(MIN_TM, 0)
        self.NUM_PRIMERS = max(NUM_PRIMERS, 0)
        self.MIN_LENGTH = max(MIN_LENGTH, 0)
        self.MAX_LENGTH = max(MAX_LENGTH, 0)
        self.COL_SIZE = max(COL_SIZE, 0)
        self.WARN_CUTOFF = max(WARN_CUTOFF, 0)

    def __repr__(self):
        return repr(self.__dict__)

    def __str__(self):
        return repr(self.__dict__)


    def get(self, key):
        key = key.upper()
        if hasattr(self, key):
            return getattr(self, key)
        elif key == 'PREFIX':
            return self.prefix
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def set(self, key, value):
        key = key.upper()
        if hasattr(self, key) or key == 'PREFIX':
            if key == 'PREFIX':
                self.prefix = str(value)
            elif isinstance(value, (int, float)) and value > 0:
                if key == 'MIN_TM':
                    self.MIN_TM = float(value)
                else:
                    setattr(self, key, int(value))
                    if self.MIN_LENGTH > self.MAX_LENGTH:
                        print('\033[93mWARNING\033[0m: \033[92mMIN_LENGTH\033[0m is greater than \033[92mMAX_LENGTH\033[0m.')
                    elif self.NUM_PRIMERS % 2:
                        print('\033[93mWARNING\033[0m: \033[92mNUM_PRIMERS\033[0m should be even number only.')
            else:
                raise ValueError('\033[41mERROR\033[0m: Illegal value \033[95m%s\033[0m for key \033[92m%s\033[0m for \033[94m%s.set()\033[0m.\n' % (value, key, self.__class__)) 
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def reset(self):
        self.prefix = 'primer'
        self.MIN_TM = 60.0
        self.NUM_PRIMERS = 0
        self.MIN_LENGTH = 15
        self.MAX_LENGTH = 60
        self.COL_SIZE = 142
        self.WARN_CUTOFF = 3


    def design(self, sequence, MIN_TM=None, NUM_PRIMERS=None, MIN_LENGTH=None, MAX_LENGTH=None, prefix=None):
        if MIN_TM is None: MIN_TM = self.MIN_TM
        if NUM_PRIMERS is None: NUM_PRIMERS = self.NUM_PRIMERS
        if MIN_LENGTH is None: MIN_LENGTH = self.MIN_LENGTH
        if MAX_LENGTH is None: MAX_LENGTH = self.MAX_LENGTH
        if prefix is None: prefix = self.prefix

        name = prefix
        sequence = RNA2DNA(sequence)
        N_BP = len(sequence)

        is_success = True
        primers = []
        warnings = []
        assembly = {}
        misprime_score = ['', '']

        try:
            print('Precalculataing Tm matrix ...')
            Tm_precalculated = precalculate_Tm(sequence)
            print('Precalculataing misprime score ...')
            (num_match_forward, num_match_reverse, best_match_forward, best_match_reverse, misprime_score_forward, misprime_score_reverse) = check_misprime(sequence)

            print('Doing dynamics programming calculation ...')
            (scores_start, scores_stop, scores_final, choice_start_p, choice_start_q, choice_stop_i, choice_stop_j, MAX_SCORE, N_primers) = dynamic_programming(NUM_PRIMERS, MIN_LENGTH, MAX_LENGTH, MIN_TM, N_BP, misprime_score_forward, misprime_score_reverse, Tm_precalculated)
            print('Doing backtracking ...')
            (is_success, primers, primer_set, warnings) = back_tracking(N_BP, sequence, scores_final, choice_start_p, choice_start_q, choice_stop_i, choice_stop_j, N_primers, MAX_SCORE, num_match_forward, num_match_reverse, best_match_forward, best_match_reverse, self.WARN_CUTOFF)

            if is_success:
                allow_forward_line = list(' ' * N_BP)
                allow_reverse_line = list(' ' * N_BP)
                for i in range(N_BP):
                    allow_forward_line[i] = str(int(min(num_match_forward[0, i] + 1, 9)))
                    allow_reverse_line[i] = str(int(min(num_match_reverse[0, i] + 1, 9)))

                misprime_score = [''.join(allow_forward_line).strip(), ''.join(allow_reverse_line).strip()]
                assembly = Assembly(sequence, primers, name, self.COL_SIZE)
                print('\033[92mSUCCESS\033[0m: Primerize 1D design() finished.\n')
            else:
                print('\033[41mFAIL\033[0m: \033[41mNO Solution\033[0m found under given contraints.\n')
        except:
            is_success = False
            print(traceback.format_exc())
            print('\033[41mERROR\033[0m: Primerize 1D design() encountered error.\n')

        params = {'MIN_TM': MIN_TM, 'NUM_PRIMERS': NUM_PRIMERS, 'MIN_LENGTH': MIN_LENGTH, 'MAX_LENGTH': MAX_LENGTH, 'N_BP': N_BP, 'COL_SIZE': self.COL_SIZE, 'WARN_CUTOFF': self.WARN_CUTOFF}
        data = {'misprime_score': misprime_score, 'assembly': assembly, 'warnings': warnings}
        return Design_1D(sequence, name, is_success, primer_set, params, data)




@jit(nopython=True, nogil=True, cache=False)
def dynamic_programming(NUM_PRIMERS, MIN_LENGTH, MAX_LENGTH, MIN_TM, N_BP, misprime_score_forward, misprime_score_reverse, Tm_precalculated):
    # could be zero, meaning user does not know.
    num_primer_sets = int(NUM_PRIMERS / 2)
    num_primer_sets_max = int(math.ceil(N_BP / float(MIN_LENGTH)))

    misprime_score_weight = 10.0
    MAX_SCORE = N_BP * 2 + 1
    MAX_SCORE += misprime_score_weight * max(numpy.amax(misprime_score_forward), numpy.amax(misprime_score_reverse)) * 2 * num_primer_sets_max

    scores_start = MAX_SCORE * numpy.ones((N_BP, N_BP, num_primer_sets_max))
    scores_stop = MAX_SCORE * numpy.ones((N_BP, N_BP, num_primer_sets_max))
    scores_final = MAX_SCORE * numpy.ones((N_BP, N_BP, num_primer_sets_max))

    # used for backtracking:
    choice_start_p = numpy.zeros((N_BP, N_BP, num_primer_sets_max))
    choice_start_q = numpy.zeros((N_BP, N_BP, num_primer_sets_max))
    choice_stop_i = numpy.zeros((N_BP, N_BP, num_primer_sets_max))
    choice_stop_j = numpy.zeros((N_BP, N_BP, num_primer_sets_max))

    # basic setup -- first primer
    # First set is special.
    #  |                     p
    #  ---------------------->
    #                   ||||||
    #                   <-----...
    #                   q
    #
    for p in range(MIN_LENGTH, MAX_LENGTH + 1):
        # STOP[reverse](1)
        q_min = max(1, p - MAX_LENGTH + 1)
        q_max = p

        for q in range(q_min, q_max + 1):
            if (Tm_precalculated[q - 1, p - 1] > MIN_TM):
                scores_stop[p - 1, q - 1, 0] = (q - 1) + 2 * (p - q + 1)
                scores_stop[p - 1, q - 1, 0] += misprime_score_weight * (misprime_score_forward[0, p - 1] + misprime_score_reverse[0, q - 1])

    best_min_score = MAX_SCORE
    n = 1
    while (n <= num_primer_sets_max):
        # final scoring -- let's see if we can 'close' at the end of the sequence.
        #
        #                 p
        #  --------------->
        #            ||||||
        #            <---------------------
        #            q                    N_BP
        #
        for p in range(1, N_BP + 1):
            q_min = max(1, p - MAX_LENGTH + 1)
            q_max = p

            # STOP[reverse]
            for q in range(q_min, q_max + 1):
                # previous primer ends had overlap with good Tm and were scored
                if (scores_stop[p - 1, q - 1, n - 1] < MAX_SCORE):
                    i = N_BP + 1
                    j = N_BP
                    last_primer_length = j - q + 1
                    if last_primer_length <= MAX_LENGTH and last_primer_length >= MIN_LENGTH:
                        scores_final[p - 1, q - 1, n - 1] = scores_stop[p - 1, q - 1, n - 1] + (i - p - 1)
                        scores_final[p - 1, q - 1, n - 1] += misprime_score_weight * (misprime_score_forward[0, p - 1] + misprime_score_reverse[0, q - 1])

        min_score = numpy.amin(scores_final[:, :, n - 1])
        if (min_score < best_min_score or n == 1):
            best_min_score = min_score
            best_n = n

        if (n >= num_primer_sets_max):
            break
        if (num_primer_sets > 0 and n == num_primer_sets):
            break

        # considering another primer set
        n += 1

        #
        #        p              i
        #  ------>              ------ ... ->
        #    |||||              ||||||
        #    <------------------------
        #    q                       j
        #
        for p in range(1, N_BP + 1):
            # STOP[forward](1)
            q_min = max(1, p - MAX_LENGTH + 1)
            q_max = p

            # STOP[reverse](1)
            for q in range(q_min, q_max + 1):
                # previous primer ends had overlap with good Tm and were scored
                if (scores_stop[p - 1, q - 1, n - 2] < MAX_SCORE):
                    # START[reverse](1)
                    min_j = max(p + 1, q + MIN_LENGTH - 1)
                    max_j = min(N_BP, q + MAX_LENGTH - 1)

                    for j in range(min_j, max_j + 1):
                        # start[reverse](2)
                        min_i = max(p + 1, j - MAX_LENGTH + 1)
                        max_i = j

                        for i in range(min_i, max_i + 1):
                            # at some PCR starge thiw will be an endpoint!
                            if (Tm_precalculated[i - 1, j - 1] > MIN_TM):
                                potential_score = scores_stop[p - 1, q - 1, n - 2] + (i - p - 1) + 2 * (j - i + 1)
                                if (potential_score < scores_start[i - 1, j - 1, n - 2]):
                                    scores_start[i - 1, j - 1, n - 2] = potential_score
                                    choice_start_p[i - 1, j - 1, n - 2] = p - 1
                                    choice_start_q[i - 1, j - 1, n - 2] = q - 1

        #
        #             i                     p
        #             ---------------------->
        #             ||||||           ||||||
        #  <----------------           <----- ...
        #                  j           q
        #
    
        # START[reverse](1)
        for j in range(1, N_BP + 1):
            # START[reverse](2)
            min_i = max(1, j - MAX_LENGTH + 1)
            max_i = j

            for i in range(min_i, max_i + 1):
                # could also just make this 1:N_BP, but that would wast a little time.
                if (scores_start[i - 1, j - 1, n - 2] < MAX_SCORE):
                    # STOP[reverse](1)
                    min_p = max(j + 1, i + MIN_LENGTH - 1)
                    max_p = min(N_BP, i + MAX_LENGTH - 1)

                    for p in range(min_p, max_p + 1):
                        # STOP[reverse](2)
                        min_q = max(j + 1, p - MAX_LENGTH + 1)
                        max_q = p

                        for q in range(min_q, max_q + 1):
                            if (Tm_precalculated[q - 1, p - 1] > MIN_TM):
                                potential_score = scores_start[i - 1, j - 1, n - 2] + (q - j - 1) + 2 * (p - q + 1)
                                potential_score += misprime_score_weight * (misprime_score_forward[0, p - 1] + misprime_score_reverse[0, q - 1])
                                if (potential_score < scores_stop[p - 1, q - 1, n - 1]):
                                    scores_stop[p - 1, q - 1, n - 1] = potential_score
                                    choice_stop_i[p - 1, q - 1, n - 1] = i - 1
                                    choice_stop_j[p - 1, q - 1, n - 1] = j - 1

    if (num_primer_sets > 0):
        N_primers = num_primer_sets
    else:
        N_primers = best_n
    return (scores_start, scores_stop, scores_final, choice_start_p, choice_start_q, choice_stop_i, choice_stop_j, MAX_SCORE, N_primers)


def back_tracking(N_BP, sequence, scores_final, choice_start_p, choice_start_q, choice_stop_i, choice_stop_j, N_primers, MAX_SCORE, num_match_forward, num_match_reverse, best_match_forward, best_match_reverse, WARN_CUTOFF):
    y = numpy.amin(scores_final[:, :, N_primers - 1], axis=0)
    idx = numpy.argmin(scores_final[:, :, N_primers - 1], axis=0)
    min_scroe = numpy.amin(y)
    q = numpy.argmin(y)
    p = idx[q]

    is_success = True
    primer_set = []
    misprime_warn = []
    primers = numpy.zeros((3, 2 * N_primers))
    if (min_scroe == MAX_SCORE):
        is_success = False
    else:
        primers[:, 2 * N_primers - 1] = [q, N_BP - 1, -1]
        for m in range(N_primers - 1, 0, -1):
            i = choice_stop_i[p, q, m]
            j = choice_stop_j[p, q, m]
            primers[:, 2 * m] = [i, p, 1]
            p = choice_start_p[i, j, m - 1]
            q = choice_start_q[i, j, m - 1]
            primers[:, 2 * m - 1] = [q, j, -1]
        primers[:, 0] = [0, p, 1]
        primers = primers.astype(int)

        for i in range(2 * N_primers):
            primer_seq = sequence[primers[0, i]:primers[1, i] + 1]
            if primers[2, i] == -1:
                primer_set.append(reverse_complement(primer_seq))

                # mispriming "report"
                end_pos = primers[0, i]
                if (num_match_reverse[0, end_pos] >= WARN_CUTOFF):
                    problem_primer = find_primers_affected(primers, best_match_reverse[0, end_pos])
                    misprime_warn.append((i + 1, num_match_reverse[0, end_pos] + 1, best_match_reverse[0, end_pos] + 1, problem_primer))
            else:
                primer_set.append(str(primer_seq))

                # mispriming "report"
                end_pos = primers[1, i]
                if (num_match_forward[0, end_pos] >= WARN_CUTOFF):
                    problem_primer = find_primers_affected(primers, best_match_forward[0, end_pos])
                    misprime_warn.append((i + 1, num_match_forward[0, end_pos] + 1, best_match_forward[0, end_pos] + 1, problem_primer))

    return (is_success, primers, primer_set, misprime_warn)


def find_primers_affected(primers, pos):
    primer_list = []
    for i in range(primers.shape[1]):
        if (pos >= primers[0, i] and pos <= primers[1, i]):
            primer_list.append(i + 1)
    return primer_list


def design_primers_1D(sequence, MIN_TM=None, NUM_PRIMERS=None, MIN_LENGTH=None, MAX_LENGTH=None, prefix=None):
    prm = Primerize_1D()
    res = prm.design(sequence, MIN_TM, NUM_PRIMERS, MIN_LENGTH, MAX_LENGTH, prefix)
    return res


def main():
    parser = argparse.ArgumentParser(description='\033[92mPrimerize 1D PCR Assembly Design\033[0m', epilog='\033[94mby Siqi Tian, 2016\033[0m', add_help=False)
    parser.add_argument('sequence', type=str, help='DNA Template Sequence')
    parser.add_argument('-p', metavar='prefix', type=str, help='Display Name of Construct', dest='prefix', default='primer')
    group1 = parser.add_argument_group('advanced options')
    group1.add_argument('-t', metavar='MIN_TM', type=float, help='Minimum Annealing Temperature', dest='MIN_TM', default=60.0)
    group1.add_argument('-n', metavar='NUM_PRIMERS', type=int, help='Number of Primers', dest='NUM_PRIMERS', default=0)
    group1.add_argument('-l', metavar='MIN_LENGTH', type=int, help='Minimum Length of each Primer', dest='MIN_LENGTH', default=15)
    group1.add_argument('-u', metavar='MAX_LENGTH', type=int, help='Maximum Length of each Primer', dest='MAX_LENGTH', default=60)
    group2 = parser.add_argument_group('commandline options')
    group2.add_argument('-q', '--quiet', action='store_true', dest='is_quiet', help='Suppress Results Printing to stdout')
    group2.add_argument('-f', '--file', action='store_true', dest='is_file', help='Write Results to Text File')
    group2.add_argument('-h', '--help', action='help', help='Show this Help Message')
    args = parser.parse_args()

    t0 = time.time()
    res = design_primers_1D(args.sequence, args.MIN_TM, args.NUM_PRIMERS, args.MIN_LENGTH, args.MAX_LENGTH, args.prefix)
    if res.is_success:
        if not args.is_quiet:
            print(res)
        if args.is_file:
            res.save()
    print('Time elapsed: %.1f s.' % (time.time() - t0))


if __name__ == "__main__":
    main()

