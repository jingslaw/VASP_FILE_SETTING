#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import os
import re
import numpy as np

__version__ = 2.0
__author__ = 'Weiguo Jing'


class Poscar(object):

    def __init__(self, filename=None):
        if filename:
            assert os.path.isfile(filename), '%s cannot be found!' % filename
            self.original = [line for line in open(filename) if line.strip()]
            char_flag = 0
            num_flag = 0
            for i, line in enumerate(self.original):
                if re.search('[a-zA-Z]+', line.strip()) and char_flag == 0:
                    self._cell_name = ''.join(line.split())
                    self._cell_name_marker = i
                    char_flag += 1
                elif re.search('[a-zA-Z]+', line.strip()) and char_flag == 1:
                    self._atoms_name_marker = i
                    char_flag += 1
                elif re.match('-?[0-9+]+', line.strip()) and num_flag == 0:
                    self._scaling_factor = float(line.strip())
                    self._scaling_factor_marker = i
                    num_flag += 1
                elif re.search('-?[0-9.+]+', line.strip()):
                    if num_flag == 1:
                        self._cell_base_vectors_marker = i
                        num_flag += 1
                    elif num_flag == 4:
                        self._atoms_number_marker = i
                        total = 0
                        for j in line.split():
                            total += int(j)
                        self.__total_atoms_number = total
                        num_flag += 1
                    else:
                        num_flag += 1
                elif re.search('Direct', line.strip()) or re.search('[CK]+', line.strip()):
                    self._atoms_position_marker = i + 1
                    if len(self.original) > self._atoms_position_marker + self.__total_atoms_number:
                        self._atoms_velocity_marker = self._atoms_position_marker + self.__total_atoms_number
                    break  # sometimes poscar don't have atoms velocity
        else:
            self.original = []
            self._cell_name = 'Standard\n'
            self._cell_name_marker = 0
            self._atoms_name_marker = 5
            self._scaling_factor = 1.00000000000000
            self._scaling_factor_marker = 1
            self._cell_base_vectors_marker = 2
            self._atoms_number_marker = 6
            self.__total_atoms_number = 1
            self._atoms_position_marker = 8
            self._atoms_velocity_marker = self._atoms_position_marker + self.__total_atoms_number
            self.original.append(self._cell_name)
            self.original.append('{:>19.14f}'.format(self.scaling_factor) + '     \n')
            self.original.append('     1.0000000000000000    0.0000000000000000    0.0000000000000000\n')
            self.original.append('     0.0000000000000000    1.0000000000000000    0.0000000000000000\n')
            self.original.append('     0.0000000000000000    0.0000000000000000    1.0000000000000000\n')
            self.original.append('   H\n')
            self.original.append('    1\n')
            self.original.append('Direct\n')
            self.original.append('  0.0000000000000000  0.0000000000000000  0.0000000000000000\n')
            self.original.append('  0.00000000E+00  0.00000000E+00  0.00000000E+00\n')

    @property
    def cell(self):
        if self._cell_base_vectors_marker:
            i = self._cell_base_vectors_marker
            return np.array([line.split() for line in self.original[i:i+3]], dtype=np.float64)
        else:
            return None

    @cell.setter
    def cell(self, cell):
        from numpy import require
        if not hasattr(cell, 'dtype'):
            cell = require(cell, dtype=np.float64)
        shape = list(np.shape(cell))
        if len(shape) == 1 and shape[0] == 3:
            self.original[self._cell_base_vectors_marker] = ' '
            for i in range(3):
                self.original[self._cell_base_vectors_marker] += ' {:>21.16f}'.format(cell[i])
            self.original[self._cell_base_vectors_marker] += '\n'
        elif shape[0] <= 3 and shape[1] == 3:
            for i in range(shape[0]):
                self.original[self._cell_base_vectors_marker + i] = ' '
                for j in range(3):
                    self.original[self._cell_base_vectors_marker + i] += \
                        ' {:>21.16f}'.format(cell[i][j])
                self.original[self._cell_base_vectors_marker + i] += '\n'
        else:
            print('\nERROR!Cell base vectors must be n*3 matrix!\nWhere n <= 3')

    @property
    def atoms_name_and_numbers(self):
        if self._atoms_name_marker and self._atoms_number_marker:
            d = {}
            for i in range(len(self.original[self._atoms_name_marker].split())):
                d[self.original[self._atoms_name_marker].split()[i]] =\
                    self.original[self._atoms_number_marker].split()[i]
            return d
        else:
            return None

    @atoms_name_and_numbers.setter
    def atoms_name_and_numbers(self, **kwargs):
        if len(kwargs) != 0:
            self.original[self._atoms_name_marker] = ''
            self.original[self._atoms_number_marker] = ''
            for name in kwargs:
                self.original[self._atoms_name_marker] += ' {:>5s}'.format(name)
                self.original[self._atoms_number_marker] += ' {:>5d}'.format(kwargs[name])
            self.original[self._atoms_name_marker] += '\n'
            self.original[self._atoms_number_marker] += '\n'

    @property
    def total_atoms_number(self):
        if self.__total_atoms_number:
            return self.__total_atoms_number
        else:
            return 0

    @total_atoms_number.setter
    def total_atoms_number(self, number=0):
        if number:
            self.__total_atoms_number = number
        else:
            print('\nERROR: total_number should be int type and bigger than 0\n')

    @property
    def atoms_position(self):
        if self._atoms_position_marker:
            start = self._atoms_position_marker
            end = start + self.__total_atoms_number
            test = np.array([line.split() for line in self.original[start:end]], dtype=float)
            return test
        else:
            return None

    @atoms_position.setter
    def atoms_position(self, position, replaced_single_atom_index=0):
        from numpy import require
        if not hasattr(position, 'dtype'):
            position = require(position, dtype=np.float64)
        shape = list(np.shape(position))
        if len(shape) == 1 and shape[0] == 3:
            if replaced_single_atom_index <= self.__total_atoms_number:
                num = replaced_single_atom_index
                self.original[self._atoms_position_marker + num] = ''
                for i in range(3):
                    self.original[self._atoms_position_marker + num] += ' {:>19.16f}'.format(position[i])
                self.original[self._atoms_position_marker + num] += '\n'
            else:
                print('\nERROR! the number of replaced single atom is out of range!\n')
        elif shape[0] <= self.__total_atoms_number and shape[1] == 3:
            for i in range(shape[0]):
                self.original[self._atoms_position_marker + i] = ''
                for j in range(3):
                    self.original[self._atoms_position_marker + i] += ' {:>19.16f}'.format(position[i][j])
                self.original[self._atoms_position_marker + i] += '\n'
        elif shape[1] == 3:
            temp = self.original[self._atoms_velocity_marker:len(self.original)]
            while len(self.original) > self._atoms_position_marker:
                self.original.pop()
            for i in range(shape[0]):
                new_atom_position = ''
                for j in range(3):
                    new_atom_position += ' {:>19.16f}'.format(position[i][j])
                new_atom_position += '\n'
                self.original.append(new_atom_position)
            self.original.extend(temp)
        else:
            print('\nERROR! one atom\'s position need three parameters\n')

    @property
    def atoms_velocity(self):
        if self._atoms_velocity_marker:
            start = self._atoms_velocity_marker
            end = start + self.__total_atoms_number
            return np.array([line.split() for line in self.original[start:end]], dtype=float)
        else:
            return None

    @atoms_velocity.setter
    def atoms_velocity(self, velocity, replaced_single_atom_index=0):
        from numpy import require
        if not hasattr(velocity, 'dtype'):
            velocity = require(velocity, dtype=np.float64)
        shape = list(np.shape(velocity))
        if len(shape) == 1 and shape[0] == 3:
            if replaced_single_atom_index <= self.__total_atoms_number:
                self.original[self._atoms_velocity_marker] = ''
                for i in range(3):
                    self.original[self._atoms_velocity_marker] += ' {:>15.8E}'.format(velocity[i])
                    print(self.original[self._atoms_velocity_marker])
                self.original[self._atoms_velocity_marker] += '\n'
            else:
                print('\nERROR! the number of replaced single atom is out of range!\n')
        elif shape[0] <= self.__total_atoms_number and shape[1] == 3:
            for i in range(shape[0]):
                self.original[self._atoms_velocity_marker + i] = ''
                for j in range(3):
                    self.original[self._atoms_velocity_marker + i] += ' {:>15.8e}'.format(velocity[i][j])
                self.original[self._atoms_velocity_marker + i] += '\n'
        elif shape[1] == 3:
            while len(self.original) > self._atoms_velocity_marker:
                self.original.pop()
            for i in range(shape[0]):
                new_atom_velocity = ''
                for j in range(3):
                    new_atom_velocity += ' {:>15.8e}'.format(velocity[i][j])
                new_atom_velocity += '\n'
                self.original.append(new_atom_velocity)

    @property
    def system_name(self):
        return self._cell_name

    @system_name.setter
    def system_name(self, new_name):
        self._cell_name = new_name
        self.original[self._cell_name_marker] = re.sub('[\w]+', new_name, self.original[self._cell_name_marker])

    @property
    def scaling_factor(self):
        return self._scaling_factor

    @scaling_factor.setter
    def scaling_factor(self, scaling_factor):
        self._scaling_factor = scaling_factor
        self.original[self._scaling_factor_marker] = \
            re.sub('[0-9.]+', '{:>.14f}'.format(self._scaling_factor), self.original[self._scaling_factor_marker])

    def print_poscar_as_file(self, filename='POSCAR'):
        with open(filename, 'w+') as fp:
            fp.writelines(self.original[:self._atoms_velocity_marker])
            fp.write('\n')
            fp.writelines(self.original
                          [self._atoms_velocity_marker:self._atoms_velocity_marker + self.__total_atoms_number])
            print('\nThe file has been written successfully\n')
            fp.close()
            return 1

    def copy(self):
        from copy import deepcopy
        return deepcopy(self)
