#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

from kltisystems import k1OrderLTIsysSisoContinuous, k1OrderLTIsysSisoDiscrete
from kltisystems import fn_example
from math        import exp
import numpy as np
#import math
#import pytest
#from unittest.mock import patch
#from numpy         import pi, dot
#from math          import sqrt

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class TestClass_Example:

    def test_example(self):
        fn_example()

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
def fn_solution(y0, x, pole, t):
    return x + ((y0-x) * exp(pole * t))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class TestClass_1OrderContinuos:

    def test_decay(self):
        pole = -2.0
        y0   =  3.0

        f1c = k1OrderLTIsysSisoContinuous(pole, y0)

        for t in [0, 0.5, 1, 2]:
            assert abs(f1c.update(t, 1.0) - fn_solution(y0, 1.0, pole, t)) < 1e-7

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class TestClass_1OrderDiscrete:

    def test_decay(self):
        pole = -2.0
        y0   = 3.0
        Ts   = 1./1e6 # [s]

        f1d = k1OrderLTIsysSisoDiscrete(pole, Ts, y0)

        for t in np.arange(Ts, 2.0, Ts):
            y = f1d.update(1.0)
            #print(f"y = {y};    sol = {fn_solution(y0, 1.0, pole, t)}")
            assert abs(y - fn_solution(y0, 1.0, pole, t)) < 1e-6

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
