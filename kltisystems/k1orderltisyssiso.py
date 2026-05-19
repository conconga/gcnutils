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

#################################
def fn_example_siso():

    import matplotlib.pyplot      as plt

    a    = -10.  # pole for the transfer function:
    tmax = 2.0
    Ts   = 5e-3 # sample rate
    T    = [i*Ts for i in range(int(tmax/Ts)+1)] # time vector
    f1c  = k1OrderLTIsysSisoContinuous(a, 0)
    f1d  = k1OrderLTIsysSisoDiscrete(a, Ts, 1)

    log_x  = list()
    log_yc = list()
    log_yd = list()
    for t in T:
        if (not (t%0.5)):
            x = float(np.random.rand()) # plant input

        f1c.update(t, x)
        f1d.update(x)
        log_x.append(x)
        log_yc.append(f1c.y)
        log_yd.append(f1d.y)
    

    #.............................................#
    #---- new figure:
    fig = 0

    #---- new figure:
    fig = fig + 1

    f = plt.figure(fig)
    f.canvas.draw() 
    f.canvas.flush_events()

    f, ax = plt.subplots(1,1,num=fig)

    ax.plot(T, log_x, T, log_yc, T, log_yd)
    ax.grid(True)
    ax.legend(('reference', 'continuous', 'discrete'))

    #.............................................#
    plt.show(block=False)
    #.............................................#

#################################
if __name__ == "__main__":
    fn_example_siso()
