import math
import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as pyplot
import numpy
import os
import re
import xlwt

from .thermo import calc_Tm


class Assembly(object):
    """Collection of result data essential for drawing an assembly scheme.

    Args:
        sequence: ``str``: Sequence of assembly design.
        primers: ``list(list(int)``: Numeric representation (end numbering and direction) of primers.
        name: ``str``: Construct prefix/name.
        COL_SIZE: ``int``: `(Optional)` Column width for assembly output. Positive number only.

    Attributes:
        sequence: ``str``: Sequence of assembly design.
        primers: ``list(list(int))``: Numeric representation (end numbering and direction) of primers.
        name: ``str``: Construct prefix/name.
        bp_lines: ``list(str)``: Strings for base-pairing lines (``'|'``).
        seq_lines: ``list(str)``: Strings for primer sequence lines.
        print_lines: ``list(tuple(str, str))``: Strings for all lines assembled, i.e. ``list(tuple('marker', 'print_line'))``.
        Tm_overlaps: ``list(float)``: List of melting temperature for all overlapping regions.

    """

    def __init__(self, sequence, primers, name, COL_SIZE=142):
        self.sequence = sequence
        self.primers = primers
        self.name = name
        (self.bp_lines, self.seq_lines, self.print_lines, self.Tm_overlaps) = _draw_assembly(self.sequence, self.primers, COL_SIZE)

        self.primer_set = []
        for i in xrange(self.primers.shape[1]):
            primer_seq = self.sequence[self.primers[0, i]:self.primers[1, i] + 1]
            if self.primers[2, i] == -1:
                self.primer_set.append(reverse_complement(primer_seq))
            else:
                self.primer_set.append(str(primer_seq))


    def __repr__(self):
        """Representation of the ``Assembly`` class.
        """

        return '\033[94m%s\033[0m {\n    \033[93m\'primers\'\033[0m: %s, \n    \033[93m\'seq_lines\'\033[0m: \033[91mlist\033[0m(\033[91mstring\033[0m * %d), \n    \033[93m\'bp_lines\'\033[0m: \033[91mlist\033[0m(\033[91mstring\033[0m * %d), \n    \033[93m\'print_lines\'\033[0m: \033[91mlist\033[0m(\033[91mtuple\033[0m * %d), \n    \033[93m\'Tm_overlaps\'\033[0m: %s\n}' % (self.__class__, repr(self.primers), len(self.seq_lines), len(self.bp_lines), len(self.print_lines), repr(self.Tm_overlaps))

    def __str__(self):
        """Results of the ``Assembly`` class. Calls ``echo()``.
        """

        return self.echo()


    def echo(self):
        """Print result in rich-text.

        Returns:
            ``str``
        """

        output = ''
        x = 0
        for i, (flag, string) in enumerate(self.print_lines):
            if (flag == '$' and 'xx' in string):
                Tm = '%2.1f' % self.Tm_overlaps[x]
                output += string.replace('x' * len(Tm), '\033[41m%s\033[0m' % Tm) + '\n'
                x += 1
            elif (flag == '^' or flag == '!'):
                num = ''.join(re.findall("[0-9]+", string))
                string = string.replace(num, '\033[100m%s\033[0m' % num) + '\n'
                if flag == '^':
                    string = string.replace('A', '\033[94mA\033[0m').replace('G', '\033[94mG\033[0m').replace('C', '\033[94mC\033[0m').replace('T', '\033[94mT\033[0m')
                else:
                    string = string.replace('A', '\033[95mA\033[0m').replace('G', '\033[95mG\033[0m').replace('C', '\033[95mC\033[0m').replace('T', '\033[95mT\033[0m')
                output += string
            elif (flag == '~'):
                output += '\033[92m%s\033[0m' % string + '\n'
            elif (flag == '='):
                output += '\033[96m%s\033[0m' % string + '\n'
            else:
                output += string + '\n'

        return output[:-1]


    def save(self, path='./', name=None):
        """Save result to text file.

        Args:
            path: ``str``: `(Optional)` Path for file saving. Use either relative or absolute path.
            name: ``str``: `(Optional)` Prefix/name for file name. When nonspecified, current object's name is used.
        """

        name = self.name if name is None else name
        f = open(os.path.join(path, '%s_assembly.txt' % name), 'w')
        lines = self.echo()
        lines += '\n%s%s\tSEQUENCE\n' % ('PRIMERS'.ljust(20), 'LENGTH'.ljust(10))
        for i, primer in enumerate(self.primer_set):
            name = '%s-\033[100m%s\033[0m%s' % (self.name, i + 1, _primer_suffix(i))
            lines += '%s\033[93m%s\033[0m\t%s\n' % (name.ljust(39), str(len(primer)).ljust(10), _primer_suffix(i).replace(' R', primer).replace(' F', primer))

        lines = lines.replace('\033[0m', '').replace('\033[100m', '').replace('\033[92m', '').replace('\033[93m', '').replace('\033[94m', '').replace('\033[95m', '').replace('\033[96m', '').replace('\033[41m', '')
        f.write(lines)
        f.close()



class Plate_96Well(object):
    """Abstraction of 96-well plates.

    Args:
        tag: ``int``: `(Optional)` Mutation library tag. Use **which_lib** number.

    Attributes:
        coords: ``set(str)``: Filled 96-Well Coordinates.
        _data: Data of primers and names, in format of ``dict: { str: tuple(primerize.util.Mutation, str) }``, i.e. ``dict: {'coord': ('tag', 'primer') }``.
    """

    def __init__(self, tag=1):
        self.coords = set()
        self._data = {}
        self.tag = 'Lib%d-' % tag

    def __repr__(self):
        """Representation of the ``Plate_96Well`` class.
        """

        if len(self):
            return '\033[94m%s\033[0m {\033[93m\'coords\'\033[0m: %s, \033[93m\'data\'\033[0m: \033[91mdict\033[0m(\033[91mtuple\033[0m * %d)}' % (self.__class__, ' '.join(sorted(self.coords)), len(self._data))
        else:
            return '\033[94m%s\033[0m (empty)' % self.__class__

    def __str__(self):
        """Results of the ``Plate_96Well`` class. Calls ``echo()``.
        """

        return self.echo()

    def __len__(self):
        """Number of filled wells.

        Returns:
            ``int``
        """

        return len(self.coords)

    def __contains__(self, coord):
        """Test if data of a given WellPosition is present.

        Args:
            coord: ``str``: WellPosition (e.g. ``'A01'``) for data.

        Returns:
            ``bool``
        """

        return coord in self.coords


    def get(self, coord):
        """Get data of a particular well or number of wells filled.

        Args:
            coord: ``str``: Keyword of parameter. Use WellPosition for well data.

        Returns:
            value of specified **coord**.

        Raises:
            AttributeError: For illegal WellPosition (out of ``range(0, 96) + 1``).
            KeyError: For nonexisted **coord**.
        """

        coord = _format_coord(coord)
        if coord_to_num(coord) is None:
            raise AttributeError('\033[41mERROR\033[0m: Illegal coordinate value \033[95m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (coord, self.__class__))
        elif coord in self:
            return self._data[coord_to_num(coord)]
        else:
            raise KeyError('\033[41mERROR\033[0m: Non-Existent coordinate value \033[95m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (coord, self.__class__))


    def set(self, coord, tag, primer):
        """Record data of a particular well.

        Args:
            coord: ``str``: WellPosition for data. Use same range as ``get()``. Existing data for the same well is overwritten.
            tag: ``primerize.util.Mutation``: Mutant representd by ``primerize.util.Mutation``. ``str`` is only supported for backward compatibility.
            primer: ``str``: Primer seuqence of well. Use sense-strand.

        Raises:
            AttributeError: For illegal WellPosition.
        """

        coord = _format_coord(coord)
        if coord_to_num(coord) is None:
            raise AttributeError('\033[41mERROR\033[0m: Illegal coordinate value \033[95m%s\033[0m for \033[94m%s.set()\033[0m.\n' % (coord, self.__class__))
        else:
            self.coords.add(coord)
            self._data[coord_to_num(coord)] = (tag, primer)


    def reset(self):
        """Clear current plate data.
        """

        self.coords = set()
        self._data = {}


    def echo(self, ref_primer=''):
        """Print result in rich-text.

        Args:
            ref_primer: ``list(str)``: `(Optional)` List of Wild-type **primer_set** for highlighting. If nonspecified, highlighting is disabled.

        Returns:
            ``str``
        """

        return _print_primer_plate(self, ref_primer)


    def save(self, ref_primer='', file_name='./plate.svg', title=''):
        """Save plate layout to image file (`SVG`).

        Args:
            ref_primer: ``list(str)``: `(Optional)` List of Wild-type primer_set for highlighting. If nonspecified, highlighting is disabled.
            file_name: ``str``: `(Optional)` File name. Include path into **file_name** when specifying. Use either relative or absolute path.
            title: ``str``: `(Optional)` Title to display on image. LaTex NOT supported.
        """

        fig = pyplot.figure()
        ax = pyplot.subplot(111)
        ax.set_aspect('equal')
        pyplot.axis([0, 13.875, 0, 9.375])
        pyplot.xticks([x * 1.125 + 0.75 for x in xrange(12)], [str(x + 1) for x in xrange(12)], fontsize=14)
        pyplot.yticks([y * 1.125 + 0.75 for y in xrange(8)], list('ABCDEFGH'), fontsize=14)
        fig.suptitle(title, fontsize=16, fontweight='bold')

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
                    tag = self._data[num][0]
                    if (isinstance(tag, Mutation) and not tag) or (isinstance(tag, str) and 'WT' in tag):
                        x_green.append(j * 1.125 + 0.75)
                        y_green.append(i * 1.125 + 0.75)
                    else:
                        if ref_primer and ref_primer == self._data[num][1]:
                            x_green.append(j * 1.125 + 0.75)
                            y_green.append(i * 1.125 + 0.75)
                        else:
                            x_violet.append(j * 1.125 + 0.75)
                            y_violet.append(i * 1.125 + 0.75)
                else:
                    x_gray.append(j * 1.125 + 0.75)
                    y_gray.append(i * 1.125 + 0.75)
        ax.scatter(x_gray, y_gray, 961, c='#ffffff', edgecolor='#333333', linewidth=5)
        ax.scatter(x_violet, y_violet, 961, c='#ecddf4', edgecolor='#c28fdd', linewidth=5)
        ax.scatter(x_green, y_green, 961, c='#beebde', edgecolor='#29be92', linewidth=5)

        matplotlib.rcParams['svg.fonttype'] = 'none'
        matplotlib.rcParams['xtick.labelsize'] = 14
        matplotlib.rcParams['ytick.labelsize'] = 14
        pyplot.savefig(file_name, orientation='landscape', format='svg')
        pyplot.close(fig)



class Mutation(object):
    """Collection of mutations for a construct.

    Args:
        mut_list: ``list(str)``: `(Optional)` List of mutations. When nonspecified, an empty instance is created; when specified, it calls ``push()``. An empty instance means no mutations, i.e. Wild-type.

    Attributes:
        _data: Data of mutations, in format of ``dict: { int: tuple(str, str) }``, i.e. ``dict: {'seqpos': ('wt_char', 'mut_char') }``.
    """

    def __init__(self, mut_str=[]):
        self._data = {}
        if mut_str: self.push(mut_str)

    def __repr__(self):
        """Representation of the ``Mutation`` class.
        """

        return '\033[94m%s\033[0m' % self.__class__

    def __str__(self):
        """Results of the ``Mutation`` class. Calls ``echo()``.
        """

        return self.echo()

    def __len__(self):
        """Number of filled wells.

        Returns:
            ``int``
        """

        return len(self._data)

    def __eq__(self, other):
        """Comparison method for whether two ``Mutation`` objects contain the same set of mutations.
        """

        if isinstance(other, (str, unicode)) and other == 'WT': return len(self) == 0
        if isinstance(other, Mutation): other = other.list()
        return other in self and len(self) == len(other)

    def __iter__(self):
        """Iterator through all mutations.
        """

        for k in self._data.keys():
            yield k, self._data[k]

    def __contains__(self, mut_str):
        """Test if a list of given mutation is present.

        Args:
            mut_list: ``list(str)``: Mutations in format of ``'wt_char'``, ``'seq_pos'``, ``'mut_char'``, (e.g. ``['G13C', 'A15T']``).

        Returns:
            ``bool``
        """

        if isinstance(mut_str, (str, unicode)): mut_str = [mut_str]
        if not (mut_str or self._data): return True
        flag = False

        for mut in mut_str:
            seq_org = RNA2DNA(mut[0])
            seq_mut = RNA2DNA(mut[-1])
            seq_pos = int(mut[1:-1])

            flag = seq_pos in self._data and self._data[seq_pos] == (seq_org, seq_mut)
            if not flag: return flag

        return flag


    def push(self, mut_str):
        """Add a list of mutations.

        Args:
            mut_list: ``list(str)``: Mutations. Valid keywords are the same as ``has()``. Each ``'seq_pos'`` can only be mutated once. Conflicting mutations are overwritten and the most recent one is saved. ``'WT'`` is ignored.

        Raises:
            ValueError: For illegal **mut_str**.
        """

        if isinstance(mut_str, (str, unicode)): mut_str = [mut_str]
        for mut in mut_str:
            if mut == 'WT': continue

            seq_org = complement(complement(RNA2DNA(mut[0])))
            seq_mut = complement(complement(RNA2DNA(mut[-1])))
            seq_pos = int(mut[1:-1])
            if seq_org == seq_mut:
                raise ValueError('\033[41mERROR\033[0m: Unchanged sequence identity by \033[94mMutation\033[0m \033[92m%s\033[0m.\n' % mut)
            self._data[seq_pos] = (seq_org, seq_mut)


    def pop(self, mut_str):
        """Remove a list of mutations.

        Args:
            mut_list: ``list(str)``: Mutations. Valid keywords are the same as ``has()``. Mutations that are not present will result in a premature return with ``False``.

        Returns:
            ``bool``: Whether all mutations in **mut_list** are successfully removed.
        """

        if isinstance(mut_str, (str, unicode)): mut_str = [mut_str]
        for mut in mut_str:
            if mut in self:
                seq_pos = int(mut[1:-1])
                self._data.pop(seq_pos, None)
            else:
                return False
        return True


    def merge(self, other):
        """Merge 2 lists of mutations.

        Args:
            other: ``primerize.util.Mutation``: Another list of mutations.

        Raises:
            TypeError: For illegal **other**.
        """

        if not isinstance(other, Mutation):
            raise TypeError('\033[41mERROR\033[0m: Illegal input type for \033[94m%s.merge()\033[0m.\n' % self.__class__)
        for mut in other.list():
            self.push(mut)


    def list(self):
        """Return a list of all mutations.

        Returns:
            ``list(str)``
        """

        return map(lambda x: '%s%s%s' % (self._data[x][0], x, self._data[x][1]), sorted(self._data.keys()))


    def echo(self):
        """Print result in rich-text, delimited by ``';'``.

        Returns:
            ``str``
        """

        output = []
        for mut in self.list():
            if mut[-2:] == 'WT':
                output.append('\033[100mWT\033[0m')
            else:
                output.append('\033[96m%s\033[0m\033[93m%s\033[0m\033[91m%s\033[0m' % (mut[0], mut[1:-1], mut[-1]))

        output = ';'.join(output) if output else '\033[100mWT\033[0m'
        return output



class Construct_List(object):
    """Collection of mutant constructs. An empty isntance of ``primerize.util.Mutation`` is always initiated as the first element, i.e. a Wild-type well as the first of list.

    Attributes:
        _data: Data of constructs, in format of ``list(primerize.util.Mutation)``.
    """

    def __init__(self):
        self._data = [Mutation()]

    def __repr__(self):
        """Representation of the ``Construct_List`` class.
        """

        if len(self):
            return '\033[94m%s\033[0m {\033[91mlist\033[0m(%s * %d)}' % (self.__class__, repr(Mutation()), len(self))
        else:
            return '\033[94m%s\033[0m (empty)' % self.__class__

    def __str__(self):
        """Results of the ``Construct_List`` class. Calls ``echo()``.
        """

        return self.echo()

    def __len__(self):
        """Number of filled wells.

        Returns:
            ``int``
        """

        return len(self._data)

    def __iter__(self):
        """Iterator through all constructs.
        """

        for i in xrange(len(self._data)):
            yield self._data[i]

    def __contains__(self, mut_list):
        """Test if a list of given mutant construct is present.

        Args:
            mut_list: ``primerize.util.Mutation``: A mutant represented by ``primerize.util.Mutation``.

        Returns:
            ``bool``
        """

        if not isinstance(mut_list, Mutation): mut_list = Mutation(mut_list)
        for construct in self._data:
            if construct == mut_list: return True
        return False


    def push(self, mut_list):
        """Add a list of mutations.

        Args:
            mut_list: ``primerize.util.Mutation``: Mutations. A mutant represented by ``primerize.util.Mutation``. If the mutant is already present, it will return ``False``.

        Returns:
            ``bool``: Whether **mut_list** is successfully added.
        """

        if not isinstance(mut_list, Mutation): mut_list = Mutation(mut_list)
        if mut_list in self: return False
        self._data.append(mut_list)
        return True


    def pop(self, mut_list):
        """Remove a list of mutations.

        Args:
            mut_list: ``primerize.util.Mutation``: A mutant represented by ``primerize.util.Mutation``. Mutant that is not present will result in a premature return with ``False``.

        Returns:
            ``bool``: Whether **mut_list** is successfully removed.
        """

        if not isinstance(mut_list, Mutation): mut_list = Mutation(mut_list)
        for i, construct in enumerate(self._data):
            if construct == mut_list:
                self._data.pop(i)
                return True
        return False


    def merge(self, other):
        """Merge 2 lists of constructs.

        Args:
            other: ``primerize.util.Construct_List``: Another list of constructs.

        Returns:
            ``primerize.util.Construct_List``: A list of duplicated constructs between inputs.

        Raises:
            TypeError: For illegal **other**.
        """

        if not isinstance(other, Construct_List):
            raise TypeError('\033[41mERROR\033[0m: Illegal input type for \033[94m%s.merge()\033[0m.\n' % self.__class__)
        repeated = Construct_List()
        for mut in other:
            flag = self.push(mut)
            if not flag: repeated.push(mut)
        repeated.pop('WT')
        return repeated


    def list(self):
        """Return a list of all constructs.

        Returns:
            ``list(list(str))``
        """

        return map(lambda x: x.list(), self._data)


    def echo(self, prefix=''):
        """Print result in rich-text.

        Returns:
            ``str``
        """

        output = ''
        for construct in self._data:
            output += prefix + construct.echo() + '\n'
        return output



def DNA2RNA(sequence):
    """Convert a DNA sequence input to RNA.

    Args:
        sequence: ``str``: Input DNA sequence.

    Returns:
        ``str``: String of RNA
    """

    return sequence.upper().replace('T', 'U')


def RNA2DNA(sequence):
    """Convert a RNA sequence input to DNA.

    Args:
        sequence: ``str``: Input RNA sequence.

    Returns:
        ``str``: String of DNA
    """

    return sequence.upper().replace('U', 'T')


def complement(sequence):
    """Convert a DNA sequence input to its complement strand.

    Args:
        sequence: ``str``: Input DNA sequence.

    Returns:
        ``str``: String of complement DNA strand.

    Raises:
        ValueError: For illegal **sequence**.
    """

    rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'}
    try:
        sequence = map(lambda x: rc_dict[x], list(sequence))
    except KeyError:
        raise ValueError('\033[41mERROR\033[0m: Illegal sequence value \033[95m%s\033[0m for \033[94mcomplement()\033[0m.\n' % sequence)

    return ''.join(sequence)


def reverse(sequence):
    """Convert a DNA sequence input to its reverse order.

    Args:
        sequence: ``str``: Input DNA sequence.

    returns:
        ``str``: String of reverse DNA strand.
    """

    return sequence[::-1]


def reverse_complement(sequence):
    """Convert a DNA sequence input to its reverse complement strand.

    Args:
        sequence: ``str``: Input DNA sequence.

    Returns:
        ``str``: String of reverse complement DNA strand.
    """

    return complement(reverse(sequence))


def _primer_suffix(num):
    if num % 2:
        return '\033[95m R\033[0m'
    else:
        return '\033[94m F\033[0m'


def _draw_assembly(sequence, primers, COL_SIZE):
    N_primers = primers.shape[1]
    seq_line_prev = list(' ' * max(len(sequence), COL_SIZE))
    bp_lines = []
    seq_lines = []
    Tms = []

    for i in xrange(N_primers):
        primer = primers[:, i]
        seq_line = list(' ' * max(len(sequence), COL_SIZE))
        (seg_start, seg_end, seg_dir) = primer

        if (seg_dir == 1):
            seq_line[seg_start: seg_end + 1] = sequence[seg_start: seg_end + 1]

            if (seg_end + 1 < len(sequence)):
                seq_line[seg_end + 1] = '-'
            if (seg_end + 2 < len(sequence)):
                seq_line[seg_end + 2] = '>'

            num_txt = '%d' % (i + 1)
            if (seg_end + 4 + len(num_txt) < len(sequence)):
                offset = seg_end + 4
                seq_line[offset:(offset + len(num_txt))] = num_txt
        else:
            seq_line[seg_start: seg_end + 1] = map(reverse_complement, sequence[seg_start: seg_end + 1])

            if (seg_start - 1 >= 0):
                seq_line[seg_start - 1] = '-'
            if (seg_start - 2 >= 0):
                seq_line[seg_start - 2] = '<'

            num_txt = '%d' % (i + 1)
            if (seg_start - 3 - len(num_txt) >= 0):
                offset = seg_start - 3 - len(num_txt)
                seq_line[offset:(offset + len(num_txt))] = num_txt

        bp_line = list(' ' * max(len(sequence), COL_SIZE))
        overlap_seq = ''
        last_bp_pos = 1
        for j in xrange(len(sequence)):
            if (seq_line_prev[j] in 'ACGT' and seq_line[j] in 'ACGT'):
                bp_line[j] = '|'
                last_bp_pos = j
                overlap_seq += sequence[j]

        if (last_bp_pos > 1):
            Tm = calc_Tm(overlap_seq, 0.2e-6, 0.1, 0.0015)
            Tms.append(Tm)
            Tm_txt = '%2.1f' % Tm
            offset = last_bp_pos + 2
            bp_line[offset:(offset + len(Tm_txt))] = 'x' * len(Tm_txt)

        bp_lines.append(''.join(bp_line))
        seq_lines.append(''.join(seq_line))
        seq_line_prev = seq_line

    print_lines = []
    for i in xrange(int(math.floor((len(sequence) - 1) / COL_SIZE)) + 1):
        start_pos = COL_SIZE * i
        end_pos = min(COL_SIZE * (i + 1), len(sequence))
        out_line = sequence[start_pos:end_pos]
        print_lines.append(('~', out_line))

        for j in xrange(len(seq_lines)):
            if (len(bp_lines[j][end_pos:].replace(' ', '')) and ('|' not in bp_lines[j][end_pos:].replace(' ', '')) and (not len(bp_lines[j][:start_pos].replace(' ', '')))):
                bp_line = bp_lines[j][start_pos:].rstrip()
            elif ('|' not in bp_lines[j][start_pos:end_pos]):
                bp_line = ' ' * (end_pos - start_pos + 1)
            else:
                bp_line = bp_lines[j][start_pos:end_pos]
            seq_line = seq_lines[j][start_pos:end_pos]

            if len(bp_line.replace(' ', '')) or len(seq_line.replace(' ', '')):
                print_lines.append(('$', bp_line))
                print_lines.append(('^!'[j % 2], seq_line))
        print_lines.append(('$', ' ' * (end_pos - start_pos + 1)))
        print_lines.append(('=', complement(out_line)))
        print_lines.append(('', '\n'))

    return (bp_lines, seq_lines, print_lines, Tms)


def _format_coord(coord):
    coord = re.findall('^([A-H]{1}[0-9]{1,2})$', coord.upper().strip())
    if not coord: return None
    coord = coord[0]
    return coord[0] + coord[1:].zfill(2)


def coord_to_num(coord):
    """Convert a 96-Well Coordinate string to number (1-based).

    Args:
        coord: ``str``: Input WellPosition coordinate string, e.g. ``'A01'``.

    Returns:
        ``int`` or ``None`` if illegal input.
    """

    if not isinstance(coord, str): return None
    coord = re.findall('^([A-H]){1}(0[1-9]|1[0-2]){1}$', coord.upper().strip())
    if not coord: return None
    coord = ''.join(coord[0])
    row = 'ABCDEFGH'.find(coord[0])
    col = int(coord[1:])
    return (col - 1) * 8 + row + 1


def num_to_coord(num):
    """Convert a 96-Well Coordinate number (1-based) to string.

    Args:
        num: ``int``: Input WellPosition coordinate number, e.g. ``96``.

    Returns:
        ``str`` or ``None`` if illegal input.
    """

    if num < 1 or num > 96 or (not isinstance(num, int)): return None
    row = 'ABCDEFGH'[(num - 1) % 8]
    col = (num - 1) / 8 + 1
    return '%s%0*d' % (row, 2, col)


def get_mut_range(mut_start, mut_end, offset, sequence):
    """Validate and calculate mutation range based on input sequence and offset. If mutation range exceeds possible range, the maximum possible range is returned.

    Args:
        mut_start: ``int``: Lower limit of mutation range, should be based on **offset**.
        mut_end: ``int``: Upper limit of mutation range, should be based on **offset**.
        offset: ``int``: Index numbering offset.
        sequence: ``str``: The sequence (length used).

    Returns:
        ``(which_muts, mut_start, mut_end)``

        - **which_muts** - ``list(int)``: The final range of mutations.
        - **mut_start** - ``int``: The valid **mut_start**.
        - **mut_end** - ``int``: The valid **mut_end**.
    """

    if (not mut_start) or (mut_start is None): mut_start = 1 - offset
    mut_start = min(max(mut_start, 1 - offset), len(sequence) - offset)
    if (not mut_end) or (mut_end is None): mut_end = len(sequence) - offset
    mut_end = max(min(mut_end, len(sequence) - offset), 1 - offset)
    which_muts = list(range(mut_start, mut_end + 1))
    return (which_muts, mut_start, mut_end)


def _get_primer_index(primer_set, sequence):
    N_primers = len(primer_set)
    coverage = numpy.zeros((1, len(sequence)))
    primers = numpy.zeros((3, N_primers))

    for n in xrange(N_primers):
        primer = RNA2DNA(primer_set[n])
        if n % 2:
            i = sequence.find(reverse_complement(primer))
        else:
            i = sequence.find(primer)
        if i == -1:
            return ([], False)
        else:
            start_pos = i
            end_pos = i + len(primer_set[n]) - 1
            seq_dir = math.copysign(1, 0.5 - n % 2)
            primers[:, n] = [start_pos, end_pos, seq_dir]
            coverage[0, start_pos:(end_pos + 1)] = 1

    return (primers.astype(int), coverage.all())


def get_mutation(nt, lib):
    """Mutate a single nucleotide.

    Args:
        nt: ``str``: The nucleotide of interest.
        lib: ``int``: The mutation library choice; choose from (``1``, ``2``, ``3``, ``4``)::

            * 1 represents "A->U, U->A, C->G, G->C",
            * 2 represents "A->C, U->C, C->A, G->A",
            * 3 represents "A->G, U->G, C->U, G->U",
            * 4 represents "A->C, U->G, C->A, G->U",
            * 5 represents "A->C, U->G, C->G, G->C".

    Returns:
        ``str``

    Raises:
        ValueError: For illegal **lib** input.
    """

    libs = {1: 'TAGC', 2: 'CCAA', 3: 'GGTT', 4: 'CGAT', 5: 'CGGC'}
    if lib not in libs:
        raise ValueError('\033[41mERROR\033[0m: Illegal value \033[95m%s\033[0m for params \033[92mwhich_lib\033[0m.\n' % lib)
    else:
        idx = 'ATCG'.find(nt)
        return libs[lib][idx]


def _print_primer_plate(plate, ref_primer):
    if not plate: return '(empty)\n'
    string = ''
    for key in sorted(plate._data):
        string += '\033[94m%s\033[0m' % num_to_coord(key).ljust(5)
        mut = plate._data[key][0]
        if isinstance(mut, Mutation):
            offset = 20 if mut else 30
            lbl = mut.echo()
            string += plate.tag + lbl.ljust(max(offset + 27 * len(mut), len(lbl)))
        else:
            if mut[-2:] == 'WT':
                string += ('%s\033[100m%s\033[0m' % (mut[:-2], mut[-2:])).ljust(28)
            else:
                string += ('%s\033[96m%s\033[0m\033[93m%s\033[0m\033[91m%s\033[0m' % (mut[:5], mut[5], mut[6:-1], mut[-1])).ljust(50)

        if ref_primer:
            for i, ref in enumerate(ref_primer):
                if ref != plate._data[key][1][i]:
                    string += '\033[41m%s\033[0m' % plate._data[key][1][i]
                else:
                    string += plate._data[key][1][i]
        else:
            string += plate._data[key][1]
        string += '\n'

    return string


def _save_plate_layout(plates, ref_primer=[], prefix='', path='./'):
    for k in xrange(len(plates[0])):
        for p in xrange(len(plates)):
            primer_sequences = plates[p][k]
            num_primers_on_plate = len(primer_sequences)

            if num_primers_on_plate:
                if num_primers_on_plate == 1 and 'A01' in primer_sequences:
                    tag = primer_sequences.get('A01')[0]
                    if (isinstance(tag, Mutation) and not tag) or (isinstance(tag, str) and 'WT' in tag): continue

                file_name = os.path.join(path, '%s_plate_%d_primer_%d.svg' % (primer_sequences.tag[:-1], k + 1, p + 1))
                print('Creating plate image: \033[94m%s\033[0m.' % file_name)
                title = '%s_plate_%d_primer_%d' % (prefix, k + 1, p + 1)
                primer_sequences.save(ref_primer[p], file_name, title)


def _save_construct_key(keys, name, path='./', prefix=''):
    prefix = 'Lib%s-' % prefix if prefix else ''
    print('Creating keys file ...')
    lines = keys.echo(prefix)
    lines = lines.replace('\033[100m', '').replace('\033[96m', '').replace('\033[93m', '').replace('\033[91m', '').replace('\033[0m', '')
    open(os.path.join(path, '%s_keys.txt' % name), 'w').write(lines)


def _save_structures(structures, name, path='./'):
    print('Creating structures file ...')
    lines = '\n'.join(structures)
    open(os.path.join(path, '%s_structures.txt' % name), 'w').write(lines)


def _save_plates_excel(plates, ref_primer=[], prefix='', path='./'):
    for k in xrange(len(plates[0])):
        file_name = os.path.join(path, '%s_plate_%d.xls' % (prefix, k + 1))
        print('Creating plate file: \033[94m%s\033[0m.' % file_name)
        workbook = xlwt.Workbook()

        for p in xrange(len(plates)):
            primer_sequences = plates[p][k]
            num_primers_on_plate = len(primer_sequences)

            if num_primers_on_plate:
                if num_primers_on_plate == 1 and 'A01' in primer_sequences:
                    tag = primer_sequences.get('A01')[0]
                    if (isinstance(tag, Mutation) and not tag) or (isinstance(tag, str) and 'WT' in tag): continue

                sheet = workbook.add_sheet('primer_%d' % (p + 1))
                sheet.col(1).width = 256 * 15
                sheet.col(2).width = 256 * 75

                sheet.write(0, 0, 'WellPosition', xlwt.easyxf('font: bold 1'))
                sheet.write(0, 1, 'Name', xlwt.easyxf('font: bold 1'))
                sheet.write(0, 2, 'Sequence', xlwt.easyxf('font: bold 1'))
                sheet.write(0, 3, 'Notes', xlwt.easyxf('font: bold 1'))

                for i, row in enumerate(sorted(primer_sequences._data)):
                    tag = primer_sequences._data[row][0]
                    primer = primer_sequences._data[row][1]
                    if isinstance(tag, Mutation):
                        format = 'font: color blue,' if (not tag or primer == ref_primer[p]) else 'font: color black,'
                        tag = ';'.join(tag.list()) if tag else 'WT'
                        tag = primer_sequences.tag + tag
                    else:
                        format = 'font: color blue,' if 'WT' in tag else 'font: color black,'

                    sheet.write(i + 1, 0, num_to_coord(row), xlwt.easyxf(format + '  italic 1'))
                    sheet.write(i + 1, 1, tag, xlwt.easyxf(format))
                    sheet.write(i + 1, 2, primer, xlwt.easyxf(format))

        if len(workbook._Workbook__worksheets): workbook.save(file_name)


def str_to_bps(structure, offset=0):
    """Convert a dot-bracket secondary structure into base-pair tuples.

    Args:
        structure: ``str``: Input secondary struture.
        offset: ``int``: `(Optional)` Index numbering offset for output numbers.

    Returns:
        ``list(list(tuple(int, int)))``: List of helices, and each helix is a list of tuple of base-pairs with their ``seqpos``.
    """

    (lbs, lbs_pk, bps, bps_pk, helices) = ([], [], [], [], [])

    for i, char in enumerate(structure):
        if char == '(':
            lbs.append(i + 1 - offset)
        elif char == ')':
            bps.append((lbs[-1], i + 1 - offset))
            lbs.pop(-1)
        elif char == '[':
            lbs_pk.append(i + 1 - offset)
        elif char == ']':
            bps_pk.append((lbs_pk[-1], i + 1 - offset))
            lbs_pk.pop(-1)
        else:
            if len(bps):
                helices.append(sorted(bps, key=lambda x: x[0]))
                bps = []
            elif len(bps_pk):
                helices.append(sorted(bps_pk, key=lambda x: x[0]))
                bps_pk = []

    if lbs or lbs_pk:
        raise ValueError('\033[41mERROR\033[0m: Unbalanced \033[92mstructure\033[0m "\033[95m%s\033[0m".\n' % structure)
    return helices


def diff_bps(structures, offset=0):
    """Find base-pairs that are not present in all secondary structure inputs. Each input secondary structure is compared to all the others.

    Args:
        structures: ``list(str)``: Input secondary structures.
        offst: ``int``: `(Optional)` Index numbering offset for output numbers.

    Returns:
        ``list(list(tuple(int, int)))``: List of helices, and each helix is a list of tuple of base-pairs with their ``seqpos``.
    """

    if isinstance(structures, str): structures = [structures]

    if len(structures) == 1:
        return str_to_bps(structures[0], offset)
    else:
        helix_all = [helix for structure in structures for helix in str_to_bps(structure, offset)]
        bps_all = ['%d@%d' % (bp[0], bp[1]) for helix in helix_all for bp in helix]
        bps = filter(lambda x: (bps_all.count(x) < len(structures)), set(bps_all))
        bps = [(int(x[0]), int(x[1])) for x in map(lambda x: x.split('@'), bps)]

        helix_all = [filter(lambda x: x in bps, helix) for helix in helix_all]
        return filter(len, helix_all)


def _mutate_primers(plates, primers, primer_set, offset, constructs, which_lib=1, is_fillWT=False):
    for i, mut in enumerate(constructs):
        plate_num = int(math.floor(i / 96.0))
        plate_pos = i % 96 + 1
        well_tag = num_to_coord(plate_pos)
        # keep track of unmatched mutations
        is_valid = mut.list()

        for p in xrange(len(primer_set)):
            wt_primer = primer_set[p]
            if mut == 'WT':
                well_name = 'Lib%d-%s' % (which_lib, 'WT')
                plates[p][plate_num].set(well_tag, mut, wt_primer)
                continue

            mut_primer = reverse_complement(wt_primer) if primers[2, p] == -1 else wt_primer
            for k, seq in mut:
                k = k + offset - 1
                if (k >= primers[0, p] and k <= primers[1, p]):
                    m_shift = int(k - primers[0, p])
                    mut_primer = list(mut_primer)
                    # check for mismatch in mutations
                    if seq[0] != mut_primer[m_shift]:
                        raise ValueError('\033[41mERROR\033[0m: Mismatch of \033[94mMutation\033[0m %s with input sequence "\033[95m%s\033[0m" at positon \033[92m%d\033[0m.\n' % (mut, mut_primer[m_shift], k))
                    mut_primer[m_shift] = seq[1]
                    mut_primer = ''.join(mut_primer)

                    seq = '%s%d%s' % (seq[0], k - offset + 1, seq[1])
                    if seq in is_valid: is_valid.remove(seq)

            mut_primer = reverse_complement(mut_primer) if primers[2, p] == -1 else mut_primer
            if mut_primer != wt_primer or is_fillWT:
                well_name = 'Lib%d-%s' % (which_lib, mut) if mut_primer != wt_primer else 'Lib%d-%s' % (which_lib, 'WT')
                plates[p][plate_num].set(well_tag, mut, mut_primer)

        for m in is_valid:
            print('\033[93mWARNING\033[0m: Unmatched (or out of bound) \033[94mMutation\033[0m position \033[92m%s\033[0m in construct %s.\n' % (m, mut))

    return plates



