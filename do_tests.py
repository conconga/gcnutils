from ksosode.kSosodeIntegrator import *
from ksosode.kSosode           import *

if __name__ == "__main__":
    tests = kSosodeIntegratorTests()
    tests.do_tests()

    tests = kSosodeTests()
    tests.do_tests()
