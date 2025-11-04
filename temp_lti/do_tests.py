from kltisystems.k2orderltisyssiso  import *
from kltisystems.k2orderltisysmimo  import *
from kltisystems.kNOrderDerivativeSiso import *

if __name__ == "__main__":
    tests = k2OrderLTIsysSisoTests()
    tests.do_tests()

    tests = k2OrderLTIsysMimoTests()
    tests.do_tests()

    tests = kNOrderDerivativeSisoTests()
    tests.do_tests()
