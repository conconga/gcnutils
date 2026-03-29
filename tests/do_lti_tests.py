import sys    
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

# when this is called from within the folder 'tests/', we
# need to add the path to the package __init__.py:
sys.path.append("..")

from kltisystems import k2OrderLTIsysSisoTests
from kltisystems import k2OrderLTIsysMimoTests
from kltisystems import kNOrderDerivativeSisoTests

if __name__ == "__main__":
    tests = k2OrderLTIsysSisoTests()
    tests.do_tests()

    tests = k2OrderLTIsysMimoTests()
    tests.do_tests()

    tests = kNOrderDerivativeSisoTests()
    tests.do_tests()
