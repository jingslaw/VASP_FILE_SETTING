#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import numpy as np

__version__ = 1.1
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
                    self.cell_name = ''.join(line.split())
                    self.cell_name_marker = i
                    char_flag += 1
                elif re.search('[a-zA-Z]+', line.strip()) and char_flag == 1:
                    self.atoms_name_marker = i
                    char_flag += 1
                elif re.match('-?[0-9+]+', line.strip()) and num_flag == 0:
                    self.scaling_factor = float(line.strip())
                    self.scaling_factor_marker = i
                    num_flag += 1
                elif re.search('-?[0-9.+]+', line.strip()):
                    if num_flag == 1:
                        self.cell_base_vectors_marker = i
                        num_flag += 1
                    elif num_flag == 4:
                        self.atoms_number_marker = i
                        sum = 0
                        for j in line.split():
                            sum += int(j)
                        self.__total_atoms_number = sum
                        num_flag += 1
                    else:
                        num_flag += 1
                elif re.search('Direct', line.strip()) or re.search('[CK]+', line.strip()):
                    self.atoms_position_marker = i + 1
                    self.atoms_velocity_marker = self.atoms_position_marker + self.__total_atoms_number
                    break  # sometimes poscar don't have atoms velocity
        else:
            self.original = []
            self.cell_name = 'None\n'
            self.cell_name_marker = 0
            self.atoms_name_marker = 5
            self.scaling_factor = 1.00000000000000
            self.scaling_factor_marker = 1
            self.cell_base_vectors_marker = 2
            self.atoms_number_marker = 6
            self.__total_atoms_number = 1
            self.atoms_position_marker = 8
            self.atoms_velocity_marker = self.atoms_position_marker + self.__total_atoms_number
            self.original.append(self.cell_name)
            self.original.append('{:>19.14f}'.format(self.scaling_factor) + '     \n')
            self.original.append('     1.0000000000000000    0.0000000000000000    0.0000000000000000\n')
            self.original.append('     0.0000000000000000    1.0000000000000000    0.0000000000000000\n')
            self.original.append('     0.0000000000000000    0.0000000000000000    1.0000000000000000\n')
            self.original.append('\n')
            self.original.append('\n')
            self.original.append('Direct\n')
            self.original.append('  0.0000000000000000  0.0000000000000000  0.0000000000000000\n')
            self.original.append('  0.00000000E+00  0.00000000E+00  0.00000000E+00\n')

    def get_cell_base_vectors(self):
        if self.cell_base_vectors_marker:
            i = self.cell_base_vectors_marker
            return np.array([line.split() for line in self.original[i:i+3]], dtype=float)
        else:
            return None

    def get_atoms_name_and_number(self):
        if self.atoms_name_marker and self.atoms_number_marker:
            d = {}
            for i in range(len(self.original[self.atoms_name_marker].split())):
                d[self.original[self.atoms_name_marker].split()[i]] = self.original[self.atoms_number_marker].split()[i]
            return d
        else:
            return None

    def get_total_atoms_number(self):
        if self.__total_atoms_number:
            return self.__total_atoms_number
        else:
            return 0

    def get_atoms_position(self):
        if self.atoms_position_marker:
            start = self.atoms_position_marker
            end = start + self.__total_atoms_number
            return np.array([line.split() for line in self.original[start:end]], dtype=float)
        else:
            return None

    def get_atoms_velocity(self):
        if self.atoms_velocity_marker:
            start = self.atoms_velocity_marker
            end = start + self.__total_atoms_number
            return np.array([line.split() for line in self.original[start:end]], dtype=float)
        else:
            return None

    def set_atoms_total_number(self, number=None):
        if number:
            self.__total_atoms_number = number
        else:
            return None

    def substitute(self, *, cell_name=None, scaling_factor=None, cell_base_vectors=np.array([]),
                   atoms_name_and_numbers=None, atoms_position=np.array([]), atoms_velocity=np.array([]),
                   replaced_single_atom_number=0):
        if cell_name:
            self.cell_name = cell_name
            self.original[self.cell_name_marker] = re.sub('[\w]+', cell_name, self.original[self.cell_name_marker])
            return self
        if scaling_factor:
            self.scaling_factor = scaling_factor
            self.original[self.scaling_factor_marker] = \
                re.sub('[0-9.]+', '{:>.14f}'.format(self.scaling_factor), self.original[self.scaling_factor_marker])
            return self
        if cell_base_vectors.any():
            if len(np.shape(cell_base_vectors)) == 1 and np.shape(cell_base_vectors)[0] == 3:
                self.original[self.cell_base_vectors_marker] = ' '
                for i in range(3):
                    self.original[self.cell_base_vectors_marker] += ' {:>21.16f}'.format(cell_base_vectors[i])
                self.original[self.cell_base_vectors_marker] += '\n'
            elif np.shape(cell_base_vectors)[0] <= 3 and np.shape(cell_base_vectors)[1] == 3:
                for i in range(np.shape(cell_base_vectors)[0]):
                    self.original[self.cell_base_vectors_marker + i] = ' '
                    for j in range(3):
                        self.original[self.cell_base_vectors_marker + i] +=\
                            ' {:>21.16f}'.format(cell_base_vectors[i][j])
                    self.original[self.cell_base_vectors_marker + i] += '\n'
            else:
                print('\nERROR!Cell base vectors must be n*3 matrix!\nWhere n <= 3')
                return self
            return self
        if atoms_name_and_numbers:
            self.original[self.atoms_name_marker] = ''
            self.original[self.atoms_number_marker] = ''
            for s in atoms_name_and_numbers:
                self.original[self.atoms_name_marker] += ' {:>5s}'.format(s)
                self.original[self.atoms_number_marker] += ' {:>5d}'.format(atoms_name_and_numbers[s])
            self.original[self.atoms_name_marker] += '\n'
            self.original[self.atoms_number_marker] += '\n'
            return self
        if atoms_position.any():
            if len(np.shape(atoms_position)) == 1 and np.shape(atoms_position)[0] == 3:
                if replaced_single_atom_number <= self.__total_atoms_number:
                    num = replaced_single_atom_number
                    self.original[self.atoms_position_marker + num] = ''
                    for i in range(3):
                        self.original[self.atoms_position_marker + num] += ' {:>19.16f}'.format(atoms_position[i])
                    self.original[self.atoms_position_marker + num] += '\n'
                else:
                    print('\nERROR! the number of replaced single atom is out of range!\n')
            elif np.shape(atoms_position)[0] <= self.__total_atoms_number and np.shape(atoms_position)[1] == 3:
                for i in range(np.shape(atoms_position)[0]):
                    self.original[self.atoms_position_marker + i] = ''
                    for j in range(3):
                        self.original[self.atoms_position_marker + i] += ' {:>19.16f}'.format(atoms_position[i][j])
                    self.original[self.atoms_position_marker + i] += '\n'
            elif np.shape(atoms_position)[1] == 3:
                temp = self.original[self.atoms_velocity_marker:len(self.original)]
                while len(self.original) > self.atoms_position_marker:
                    self.original.pop()
                for i in range(np.shape(atoms_position)[0]):
                    new_atom_position = ''
                    for j in range(3):
                        new_atom_position += ' {:>19.16f}'.format(atoms_position[i][j])
                    new_atom_position += '\n'
                    self.original.append(new_atom_position)
                self.original.extend(temp)
            else:
                print('\nERROR! one atom\'s position need three parameters\n')
            return self
        if atoms_velocity.any():
            if len(np.shape(atoms_velocity)) == 1 and np.shape(atoms_velocity)[0] == 3:
                if replaced_single_atom_number <= self.__total_atoms_number:
                    self.original[self.atoms_velocity_marker] = ''
                    for i in range(3):
                        self.original[self.atoms_velocity_marker] += ' {:>15.8E}'.format(atoms_velocity[i])
                        print(self.original[self.atoms_velocity_marker])
                    self.original[self.atoms_velocity_marker] += '\n'
                else:
                    print('\nERROR! the number of replaced single atom is out of range!\n')
            elif np.shape(atoms_velocity)[0] <= self. __total_atoms_number and np.shape(atoms_velocity)[1] == 3:
                for i in range(np.shape(atoms_velocity)[0]):
                    self.original[self.atoms_velocity_marker + i] = ''
                    for j in range(3):
                        self.original[self.atoms_velocity_marker + i] += ' {:>15.8e}'.format(atoms_velocity[i][j])
                    self.original[self.atoms_velocity_marker + i] += '\n'
            elif np.shape(atoms_velocity)[1] == 3:
                while len(self.original) > self.atoms_velocity_marker:
                    self.original.pop()
                for i in range(np.shape(atoms_velocity)[0]):
                    new_atom_velocity = ''
                    for j in range(3):
                        new_atom_velocity += ' {:>15.8e}'.format(atoms_velocity[i][j])
                    new_atom_velocity += '\n'
                    self.original.append(new_atom_velocity)
            return self

    def print_poscar_as_file(self, *, filename='POSCAR', pathname=os.getcwd()):
        filename = pathname + '\\' + filename
        with open(filename, 'w+') as fp:
            fp.writelines(self.original[:self.atoms_velocity_marker])
            fp.write('\n')
            fp.writelines(self.original
                          [self.atoms_velocity_marker:self.atoms_velocity_marker + self.__total_atoms_number])
            print('\nThe file has been write successfully\n')
            return 1
