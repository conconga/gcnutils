import sys    
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

# when this is called from within the folder 'tests/', we
# need to add the path to the package __init__.py:
sys.path.append("..")

from ksosode import kSosodeIntegratorTests
from ksosode import kSosodeTests

if __name__ == "__main__":
    tests = kSosodeIntegratorTests()
    tests.do_tests()

    tests = kSosodeTests()
    tests.do_tests()
