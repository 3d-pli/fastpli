# -*- coding: utf-8 -*-
"""
misc methods for helping
"""

import glob
import os


def pip_freeze():
    """ turns pip freeze into a string """
    try:
        from pip._internal.operations import freeze
    except ImportError:
        from pip.operations import freeze
    return '\n'.join(freeze.freeze())


def version_file_name(file_name):
    """ Versions file name with .v{i}. Returns new file name with latest i """
    file_path = os.path.dirname(file_name)
    file_name = os.path.basename(file_name)
    files = glob.glob(os.path.join(file_path, file_name + '*'))

    def in_list(i, files):
        for f in files:
            if file_name + f'.v{i}' in f:
                return True
        return False

    i = 0
    while in_list(i, files):
        i += 1

    return os.path.join(file_path, file_name + f'.v{i}')


def version_path(path, name):
    """ Versions folder name with .v{i}. Returns new path with latest i """
    folders = [
        p for p in os.listdir(path) if os.path.isdir(os.path.join(path, p))
    ]

    def in_list(i, name):
        for f in folders:
            if name + f'.v{i}' in f:
                return True
        return False

    i = 0
    while in_list(i, name):
        i += 1

    return os.path.join(path, name + f'.v{i}')
