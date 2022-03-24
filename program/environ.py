import os
import sys


def assert_ccp4():
    if "CCP4" not in os.environ:
        sys.exit("CCP4 environment variable not set")
