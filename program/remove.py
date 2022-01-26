#!/usr/bin/python3

import glob
import os
import shutil


def _main():
    dirs = glob.glob("data/af/*") + glob.glob("data/ep/*")
    for dircetory in dirs:
        _remove_result(f"{dircetory}/ccp4i/ccp4i.pdb")
        _remove_result(f"{dircetory}/modelcraft/modelcraft.cif")


def _remove_result(path):
    if not os.path.exists(path):
        directory = os.path.dirname(path)
        try:
            shutil.rmtree(directory)
        except FileNotFoundError:
            pass
        else:
            print("Removed", directory)


if __name__ == "__main__":
    _main()
