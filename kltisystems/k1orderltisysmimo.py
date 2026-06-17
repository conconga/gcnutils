#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung: These classes implement the dynamical behaviour of a first order
    LTI system MIMO.
Autor: Luciano Auguto Kruk
Erstellt am: 00.00.2017
Version: 1.0.0
Lizenz: MIT License
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##
import numpy              as np


if __name__ == "__main__":
    from k1orderltisyssiso import k1OrderLTIsysSisoDiscrete
else:
    from .k1orderltisyssiso import k1OrderLTIsysSisoDiscrete

#################################
## First Order Discrete System ##
#################################
class k1OrderLTIsysMimoDiscrete:
    """
    The system is the discrete differential implementation of:

    f(t) = x + (y0 - x).exp(pole.t)

    where:
        pole :  (float) or ([N x 1]) pole < 0 for stability
        x    :  (vector) [N x 1] inputs
        y    :  (vector) [N x 1] outputs
        y0   :  (vector) [N x 1] initial condition

    """

    def __init__(self, pole, Ts, y0):
        """
        Inputs:
            pole :  (float) or ([N x 1]) pole < 0 for stability
            y0   :  (vector) [N x 1] initial condition
            Ts   :  (float)  [s] sampling period
        """

        self.N = len(y0) if isinstance(y0, (list, tuple)) else len(y0.squeeze())

        if isinstance(pole, (list, tuple, np.ndarray)):
            assert len(pole) == self.N
        else:
            pole = [float(pole)] * self.N

        self.filters = [k1OrderLTIsysSisoDiscrete(pole[i], Ts, y0[i]) for i in range(self.N)]
        self.y       = [ self.filters[i].y for i in range(self.N) ]

    def update(self, x):
        """
        Inputs:
            x    :  (vector) [N x 1] inputs
        """

        if isinstance(x, (list, tuple)):
            y = [ self.filters[i].update(x[i]) for i in range(self.N) ]
        else:
            x = x.reshape(-1)
            y = [ self.filters[i].update(x[i]) for i in range(self.N) ]

        if isinstance(x, np.ndarray):
            self.y = np.asarray(y)
        else:
            self.y = x.__class__(y)

        return self.y

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

