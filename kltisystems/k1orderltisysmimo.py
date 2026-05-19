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
        pole :  (float) pole < 0 for stability
        x    :  (vector) [N x 1] inputs
        y    :  (vector) [N x 1] outputs
        y0   :  (vector) [N x 1] initial condition

    """

    def __init__(self, pole, Ts, y0):
        """
        Inputs:
            pole :  (float) pole < 0 for stability
            y0   :  (vector) [N x 1] initial condition
            Ts   :  (float)  [s] sampling period
        """

        self.N = len(y0) if isinstance(y0, (list, tuple)) else len(y0.squeeze())
        self.filters = [k1OrderLTIsysSisoDiscrete(pole, Ts, y0[i]) for i in range(self.N)]

    def update(self, x):
        """
        Inputs:
            x    :  (vector) [N x 1] inputs
        """

        if isinstance(x, (list, tuple)):
            y = [ self.filters[i].update(x[i]) for i in range(self.N) ]
        else:
            x = x.squeeze()
            y = [ self.filters[i].update(x[i]) for i in range(self.N) ]

        if isinstance(x, np.ndarray):
            self.y = np.asarray(y)
        else:
            self.y = x.__class__(y)

        return self.y

#################################
def fn_example_mimo():

    import matplotlib.pyplot      as plt

    pole = -10. # pole for the transfer function:
    tmax = 1.0  # [s]
    Ts   = 5e-3 # sample rate
    T    = np.arange(0,tmax, Ts)
    y0   = [1,2,3]
    f1d  = k1OrderLTIsysMimoDiscrete(pole, Ts, y0)

    lst_log = [ (0, y0) ]
    for t in T:
        y = f1d.update([0,0,0])
        lst_log.append((t, y))

    #.............................................#
    #---- new figure:
    fig = 0

    #---- new figure:
    fig = fig + 1

    f = plt.figure(fig).clf()
    f, ax = plt.subplots(len(y0),1,num=fig, sharex=True)

    for idx in range(len(y0)):
        ax[idx].plot([i[0] for i in lst_log], [i[1][idx] for i in lst_log])
        ax[idx].grid(True)
        ax[idx].legend((f'y[{idx}]',))

    f.canvas.flush_events()
    f.canvas.draw() 

    #.............................................#
    plt.show(block=False)
    #.............................................#

#################################
if __name__ == "__main__":
    fn_example_mimo()
