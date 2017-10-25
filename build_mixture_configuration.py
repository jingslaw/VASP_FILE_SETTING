#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

__version__ = 2.0
__author__ = 'Weiguo Jing'

import os
from SetPoscar import Poscar
import numpy as np


def set_eta(*args, **kwargs):
    if 'mix_method' in kwargs:
        mix_method = kwargs.pop('mix_method')
        if mix_method == 'isometric':
            middle = left = right = None
            count = 0
            if 'middle' in kwargs:
                middle = kwargs.pop('middle')
                if isinstance(middle, int) is not True:
                    print("\nERROR:The value of 'middle' should be int type!\n")
                    return None
            if 'left' in kwargs:
                left = kwargs.pop('left')
                if isinstance(left, int) is not True:
                    print("\nERROR:The value of 'left' should be int type!\n")
                    return None
                count = count + 1
            if 'right' in kwargs:
                right = kwargs.pop('right')
                if isinstance(right, int) is not True:
                    print("\nERROR:The value of 'right' should be int type!\n")
                    return None
                count = count + 1
            if len(args) != 0:
                if len(args) >= 3:
                    print('\nWARNING: Too much parameters in isometric method.\n')
                    return None
                if len(args) != count:
                    print('\nERROR: Mismatch between number of parameters in args and inserting areas!\n')
                    return None
                eta = []
                if middle:
                    step = 1.0 / (1 + abs(middle))
                    eta = eta + [(i + 1) * step for i in range(middle)]
                if right:
                    if args[count - 1] < 0:
                        step = -args[count - 1]
                    else:
                        step = args[count - 1]
                    eta = eta + [1 + (i + 1) * step for i in range(right)]
                    count = count - 1
                if left:
                    if args[count - 1] > 0:
                        step = -args[count - 1]
                    else:
                        step = args[count - 1]
                    eta = eta + [(i + 1) * step for i in range(left)]
                return sorted(eta)
            elif middle:
                step = 1.0 / (1 + abs(middle))
                eta = [(i + 1) * step for i in range(middle)]
                if left:
                    eta = eta + [-(i + 1) * step for i in range(left)]
                if right:
                    eta = eta + [1 + (i + 1) * step for i in range(right)]
                return sorted(eta)
            else:
                print("\nERROR:set_eta function needs at least the value of 'middle' parameter!\n")
        elif mix_method == 'arbitrary':
            return args
        else:
            print("\nERROR:Mix method should be 'isometric' or 'arbitrary'!\n")
            return None
    else:
        print("\nERROR:The mix method is not defined!\n")
        return None


def construct_mixed_configuration(filename1='POSCAR1', filename2='POSCAR2', *args, **kwargs):
    epsilon = 0.0001
    poscar1 = Poscar(filename1)
    poscar2 = Poscar(filename2)
    position1 = poscar1.atoms_position
    position2 = poscar2.atoms_position
    if poscar1.total_atoms_number != poscar2.total_atoms_number:
        print('\nERROR! Two POSCAR files must have same total atoms number!\n')
        return 0
    delta = position2 - position1
    delta = np.where(np.abs(delta - 1) < epsilon, delta - 1, delta)
    delta = np.where(np.abs(delta + 1) < epsilon, delta + 1, delta)
    eta = set_eta(*args, **kwargs)
    new_position = np.zeros((len(eta), poscar1.total_atoms_number, 3))
    poscars = []
    for i, j in enumerate(eta):
        new_position[i] = position1 + j * delta
        poscar1.atoms_position = new_position[i]
        poscars.append(poscar1.copy())
    return eta, poscars


def setup_as_package(filename1='POSCAR1', filename2='POSCAR2', *args, **kwargs):
    from shutil import copy
    if 'output_path' in kwargs:
        output_path = kwargs.pop('output_path')
    else:
        output_path = os.getcwd()
    if 'output_folder' in kwargs:
        output_folder = kwargs.pop('output_folder')
    else:
        output_folder = None
    other_vasp_files = []
    if 'KPOINTS' in kwargs:
        kpoints_path = kwargs.pop('KPOINTS')
        other_vasp_files.append(kpoints_path)
    if 'INCAR' in kwargs:
        incar_path = kwargs.pop('INCAR')
        other_vasp_files.append(incar_path)
    if 'POTCAR' in kwargs:
        potcar_path = kwargs.pop('POTCAR')
        other_vasp_files.append(potcar_path)
    if output_folder:
        output_path = os.path.join(output_path, output_folder)
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
    eta, poscars = construct_mixed_configuration(filename1, filename2, *args, **kwargs)
    for i, j in enumerate(eta):
        output_path = os.path.join(output_path, 'eta=' + '{:.3f}'.format(j))
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        filename = os.path.join(output_path, 'POSCAR')
        poscars[i].print_poscar_as_file(filename)
        if len(other_vasp_files) != 0:
            for path in other_vasp_files:
                assert os.path.isfile(path), '%s cannot be found!' % path
                copy(path, output_path)
        output_path = os.path.split(output_path)[0]
    return 0


if __name__ == "__main__":
    file1 = r'/home/jingslaw/PycharmProjects/BaF2Eu@4f/POSCAR'
    file2 = r'/home/jingslaw/PycharmProjects/BaF2Eu@cb/POSCAR'
    area = [-0.1, 0.1]
    order = {'left': 4, 'middle': 9, 'right': 4, 'mix_method': 'isometric', 'output_folder': 'test',
             'KPOINTS': r'/home/jingslaw/PycharmProjects/BaF2Eu@4f/KPOINTS',
             'INCAR': r'/home/jingslaw/PycharmProjects/BaF2Eu@4f/INCAR',
             'POTCAR': r'/home/jingslaw/PycharmProjects/BaF2Eu@4f/POTCAR'}
    result = setup_as_package(file1, file2, *area, **order)
    print(result)
