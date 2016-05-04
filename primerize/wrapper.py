import math
import os

if __package__ is None or not __package__:
    from util import *
else:
    from .util import *


class Design_Single(object):
    def __init__(self, init_dict):
        for key in init_dict:
            if key not in ['sequence', 'name', 'is_success', 'primer_set', 'params', 'data']:
                raise ValueError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s\033[0m.\n' % (key, self.__class__))
            key_rename = '_' + key if key in ['params', 'data'] else key
            setattr(self, key_rename, init_dict[key])

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
            name = self.name if name is None else name
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
                    p_1 = '\033[100m%d\033[0m%s' % (warning[0], _primer_suffix(warning[0] - 1))
                    p_2 = ', '.join('\033[100m%d\033[0m%s' % (x, _primer_suffix(x - 1)) for x in warning[3])
                    output += '\033[93mWARNING\033[0m: Primer %s can misprime with %d-residue overlap to position %s, which is covered by primers: %s\n' % (p_1.rjust(4), warning[1], str(int(warning[2])).rjust(3), p_2)
                return output[:-1]

            elif key == 'primer':
                output = '%s%s\tSEQUENCE\n' % ('PRIMERS'.ljust(20), 'LENGTH'.ljust(10))
                for i in range(len(self.primer_set)):
                    name = '%s-\033[100m%s\033[0m%s' % (self.name, i + 1, _primer_suffix(i))
                    output += '%s\033[93m%s\033[0m\t%s\n' % (name.ljust(39), str(len(self.primer_set[i])).ljust(10), _primer_suffix(i).replace(' R', self.primer_set[i]).replace(' F', self.primer_set[i]))
                return output[:-1]

            elif key == 'assembly':
                return self._data['assembly'].echo()
            elif not key:
                return self.echo('misprime') + '\n' + self.echo('assembly') + '\n' + self.echo('primer') + '\n\n' + self.echo('warning') + '\n'

            else:
                raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.echo()\033[0m.\n' % (key, self.__class__))
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of key \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (key, self.__class__))



class Design_Plate(object):
    def __init__(self, init_dict):
        for key in init_dict:
            if key not in ['sequence', 'name', 'is_success', 'primer_set', 'structures', 'params', 'data']:
                raise ValueError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s\033[0m.\n' % (key, self.__class__))
            key_rename = '_' + key if key in ['params', 'data'] else key
            setattr(self, key_rename, init_dict[key])
        if self.get('TYPE') == 'Mutate-and-Map':
            self._data['illustration'] = _draw_region(self.sequence, self._params)
        elif self.get('TYPE') == 'Mutation/Rescue':
            self._data['illustration'] = _draw_str_region(self.sequence, self.structures, self._data['bps'], self._params)

    def __repr__(self):
        structures = '\033[95mstructures\033[0m = %s, \n' % repr(self.structures) if self.get('TYPE') == 'Mutation/Rescue' else ''
        return '\033[94m%s\033[0m {\n\033[95msequence\033[0m = \'%s\', \n\033[95mname\033[0m = \'%s\', \n\033[95mis_success\033[0m = \033[41m%s\033[0m, \n\033[95mprimer_set\033[0m = %s, \n%s\033[95mparams\033[0m = %s, \n\033[95mdata\033[0m = {\n    \033[92m\'constructs\'\033[0m: %s, \n    \033[92m\'assembly\'\033[0m: %s, \n    \033[92m\'plates\'\033[0m: %s\n}' % (self.__class__, self.sequence, self.name, self.is_success, repr(self.primer_set), structures, repr(self._params), repr(self._data['constructs']), repr(self._data['assembly']), repr(self._data['plates']))

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
            return self._data['constructs']
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def save(self, key='', path='./', name=None):
        if self.is_success:
            if name is None: name = self.name
            key = key.lower()
            if key == 'table':
                _save_plates_excel(self._data['plates'], self.primer_set, name, path)
            elif key == 'image':
                _save_plate_layout(self._data['plates'], self.primer_set, name, path)
            elif key == 'construct':
                _save_construct_key(self._data['constructs'], name, path, self._params['which_lib'])
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
            elif key == 'region':
                return '\n'.join(self._data['illustration']['lines'])

            elif not key:
                return self.echo('assembly') + '\n\n' + self.echo('plate') + '\n\n' + self.echo('region')
            else:
                raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.echo()\033[0m.\n' % (key, self.__class__))
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of key \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (key, self.__class__))


