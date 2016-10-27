from .util import str_to_bps


def _draw_common(fragments, labels):
    (illustration_1, illustration_2, illustration_3) = ('', '', '')

    if len(fragments[0]) >= len(labels[0]):
        illustration_1 += '\033[91m' + fragments[0][0] + '\033[0m\033[40m' + fragments[0][1:] + '\033[0m'
        illustration_2 += '\033[91m|%s\033[0m' % (' ' * (len(fragments[0]) - 1))
        illustration_3 += '\033[91m%s%s\033[0m' % (labels[0], ' ' * (len(fragments[0]) - len(labels[0])))
    elif fragments[0]:
        illustration_1 += '\033[91m' + fragments[0][0] + '\033[0m\033[40m' + fragments[0][1:] + '\033[0m'
        illustration_2 += '\033[91m%s\033[0m' % (' ' * len(fragments[0]))
        illustration_3 += '\033[91m%s\033[0m' % (' ' * len(fragments[0]))

    if len(fragments[1]) >= len(labels[1]) + len(labels[2]):
        illustration_1 += '\033[44m' + fragments[1][0] + '\033[0m\033[46m' + fragments[1][1:-1] + '\033[0m\033[44m' + fragments[1][-1] + '\033[0m'
        illustration_2 += '\033[92m|%s|\033[0m' % (' ' * (len(fragments[1]) - 2))
        illustration_3 += '\033[92m%s%s%s\033[0m' % (labels[1], ' ' * (len(fragments[1]) - len(labels[1]) - len(labels[2])), labels[2])
    elif fragments[1]:
        if len(fragments[1]) >= len(labels[1]):
            illustration_1 += '\033[44m' + fragments[1][0] + '\033[0m\033[46m' + fragments[1][1:] + '\033[0m'
            illustration_2 += '\033[92m|%s\033[0m' % (' ' * (len(fragments[1]) - 1))
            illustration_3 += '\033[92m%s%s\033[0m' % (labels[1], ' ' * (len(fragments[1]) - len(labels[1])))
        else:
            illustration_1 += '\033[46m' + fragments[1] + '\033[0m'
            illustration_2 += '\033[92m%s\033[0m' % (' ' * len(fragments[1]))
            illustration_3 += '\033[92m%s\033[0m' % (' ' * len(fragments[1]))

    if len(fragments[2]) >= len(labels[3]):
        illustration_1 += '\033[40m' + fragments[2][:-1] + '\033[0m\033[91m' + fragments[2][-1] + '\033[0m'
        illustration_2 += '\033[91m%s|\033[0m' % (' ' * (len(fragments[2]) - 1))
        illustration_3 += '\033[91m%s%s\033[0m' % (' ' * (len(fragments[2]) - len(labels[3])), labels[3])
    elif fragments[2]:
        illustration_1 += '\033[40m' + fragments[2][:-1] + '\033[0m\033[91m' + fragments[2][-1] + '\033[0m'
        illustration_2 += '\033[91m%s\033[0m' % (' ' * len(fragments[2]))
        illustration_3 += '\033[91m%s\033[0m' % (' ' * len(fragments[2]))

    return (illustration_1, illustration_2, illustration_3)


def _draw_region(sequence, params):
    offset = params['offset']
    start = params['which_muts'][0] + offset - 1
    end = params['which_muts'][-1] + offset - 1
    fragments = []

    if start <= 20:
        fragments.append(sequence[:start])
    else:
        fragments.append(sequence[:10] + '......' + sequence[start - 10:start])
    if end - start <= 40:
        fragments.append(sequence[start:end + 1])
    else:
        fragments.append(sequence[start:start + 20] + '......' + sequence[end - 19:end + 1])
    if len(sequence) - end <= 20:
        fragments.append(sequence[end + 1:])
    else:
        fragments.append(sequence[end + 1:end + 11] + '......' + sequence[-10:])

    labels = ['%d' % (1 - offset), '%d' % params['which_muts'][0], '%d' % params['which_muts'][-1], '%d' % (len(sequence) - offset)]
    (illustration_1, illustration_2, illustration_3) = _draw_common(fragments, labels)
    return {'labels': labels, 'fragments': fragments, 'lines': (illustration_1, illustration_2, illustration_3)}


def _draw_str_region(sequence, structures, bps, params):
    offset = params['offset']
    start = params['which_muts'][0] + offset - 1
    end = params['which_muts'][-1] + offset - 1
    fragments = []

    fragments.append(sequence[:start])
    fragments.append(sequence[start:end + 1])
    fragments.append(sequence[end + 1:])

    labels = ['%d' % (1 - offset), '%d' % params['which_muts'][0], '%d' % params['which_muts'][-1], '%d' % (len(sequence) - offset)]
    (illustration_1, illustration_2, illustration_3) = _draw_common(fragments, labels)

    illustration_str = ''
    bps = [bp for helix in bps for bp in helix]
    for structure in structures:
        this_bps = [bp for helix in str_to_bps(structure) for bp in helix]
        this_bps = filter(lambda x: (x in bps), this_bps)
        bps = filter(lambda x: (x not in this_bps), bps)
        this_bps = [nt for bp in this_bps for nt in bp]

        for i, nt in enumerate(structure):
            if i + 1 in this_bps:
                illustration_str += '\033[41m%s\033[0m' % nt
            else:
                illustration_str += nt
        illustration_str += '\n'

    return {'labels': labels, 'fragments': fragments, 'lines': (illustration_3, illustration_2, illustration_1, illustration_str)}
