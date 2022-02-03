import os
import testsetep
import testsetmr
import testsetaf
import pipelines


def _main():
    if "CCP4" not in os.environ:
        print("CCP4 environment variable not set")
        return
    testsetep.prepare()
    testsetmr.prepare()
    testsetaf.prepare()
    pipelines.run()


if __name__ == "__main__":
    _main()
