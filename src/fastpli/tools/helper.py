import glob
import os

from .. import __version__


def pip_freeze():
    try:
        from pip._internal.operations import freeze
    except ImportError:
        from pip.operations import freeze
    return "\n".join(freeze.freeze())


def version_file_name(file_name):

    file_path = os.path.dirname(file_name)
    file_name = os.path.basename(file_name)
    files = glob.glob(os.path.join(file_path, file_name + '*'))

    def in_list(i, file):
        for f in files:
            if file_name + ".v{}".format(i) in f:
                return True
        return False

    i = 0
    while in_list(i, files):
        i += 1

    return os.path.join(file_path, file_name + ".v{}".format(i))


def version_path(path, name):
    folders = [
        p for p in os.listdir(path) if os.path.isdir(os.path.join(path, p))
    ]

    def in_list(i, name):
        for f in folders:
            if name + f".v{i}" in f:
                return True
        return False

    i = 0
    while in_list(i, name):
        i += 1

    return os.path.join(path, name + f".v{i}")
