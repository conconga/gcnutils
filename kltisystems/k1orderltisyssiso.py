#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung: These classes implement the dynamical behaviour of a first order
    LTI system SISO.
Autor: Luciano Auguto Kruk
Erstellt am: 00.00.2017
Version: 1.0.0
Lizenz: MIT License
GitHub: 
"""
#====================================#
##WWww=--  import section: --=wwWW##
import numpy              as np
import scipy.integrate    as Int

###################################
## First Order Continuous System ##
###################################
class k1OrderLTIsysSisoContinuous:
    """
    The system is the differential implementation of:

    f(t) = x + (y0 - x).exp(pole.t)
    """

    def __init__(self, pole, y0):
        """
        T(s) = -pole/(s-pole)

        Inputs:
            pole :   system pole ( Re{pole}<0 for estability )
            y0   :   initial condition
        """

        self.pole = pole
        self.y    = y0
        self.t    = 0.

    def _dydt(self, y, t, x_):
        """
        _dydt(t) = a.y(t) - a.x(t)
        """

        x = x_
        return self.pole*(y-x)

    def update(self, t, x):
        y = Int.odeint(self._dydt, self.y, [self.t, t], (x,)) # returns y[t-1] e y[t]
        self.y = y.reshape(-1)[1]
        self.t = t

        return self.y

#################################
## First Order Discrete System ##
#################################
class k1OrderLTIsysSisoDiscrete:

    def __init__(self, a, Ts, y0):
        """
        y[t] = k*( b.x[t] - b.x[t-1] + c.y[t-1] )
        with:
          b = -a.T
          c = (2+a.T)
          k = 1/(2-a.T)
        """
        if a>0:
            raise ValueError(f'the pole "a={a}" is instable')

        k = 2.-(a*Ts)

        self.b  = -a*Ts/k
        self.c  = (2.+(a*Ts))/k
        self.y  = y0
        self.x  = 0.
        self.t  = 0.
        self.Ts = Ts

    def update(self, x):
        self.y = (self.b*x) + (self.b*self.x) + (self.c*self.y)
        self.t += self.Ts
        self.x = x

        return self.y

##########################
## k1Order Siso Factory ##
##########################
def k1OrderLTIsysSisoFactory(*args, **kargs):

    Ts = kargs.get('Ts', None)

    # when Ts is defined:
    if Ts is not None:
        if Ts > 0:
            return k1OrderLTIsysSisoDiscrete(*args, **kargs)
        else:
            raise NotImplementedError("this will NEVER be implemented!")
    else:
        # when Ts is positional:
        if len(args) >= 3:
            return k1OrderLTIsysSisoDiscrete(*args, **kargs)
        else:
            return k1OrderLTIsysSisoContinuous(*args, **kargs)

#################################
