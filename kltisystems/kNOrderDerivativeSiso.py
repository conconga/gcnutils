#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung: These classes implement a linear time-invariante (LTI) filter to
    support estimation of any derivative of the input signal. This is a SISO
    implementation.
Autor: Luciano Auguto Kruk
Erstellt am: 29.09.2025
Version: 1.0.0
Lizenz: Please keep this header with the file.
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##
import numpy            as np
import scipy.signal     as ss

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kNOrderDerivativeSiso:
    """
    This objects smoothes the input signal and calculates its derivatives.

    \parameters:                     
        order    : maximum order of the derivative to estimate; eg. y''' => order=3
        pole     : poles of the system, negative floats
        Ts       : discrete interval
    """

    def __init__(self, order, pole, Ts):
        assert order > 0
        assert pole < 0
        assert Ts > 0

        self.M    = order + 1
        self.pole = pole
        self.Ts   = Ts

#################################
## ##WWww=--  main:  --=wwWW## ##
#################################
if __name__ == "__main__":

    import matplotlib.pylab as plt

    pole = -100
    order = 4
    Fs = 400 # [Hz]
    Ts = 1./Fs
    U = list()

    # time and input signal:
    for i in range(4):
        v = np.random.randn()
        U.append( [v for i in range(Fs)] )
    U = np.hstack(U)
    T = np.asarray(range(len(U))) * Ts

    # polynomial:
    den = np.poly1d( [1, -pole] )
    for i in range(order-1):
        den *= np.poly1d( [1, -pole] )
    # den = [1 ..  a1 a0]
    den  = list(den) # indexing 'den' is non intuitive
    rden = den[::-1] # reversing

    # system:
    A = np.hstack( (np.zeros((order-1, 1)), np.eye(order-1)) )
    A = np.vstack( (A, [-i for i in rden[:order]]) )
    B = np.zeros(order)
    B[-1] = den[-1]
    B = B.reshape((-1,1))
    D = 0;

    plt.figure(1).clf()
    fig,ax = plt.subplots(order, 1, num=1, sharex=True)
    ax[0].plot(T, U)
    for i in range(order):
        C = np.zeros(order)
        C[i] = 1.0
        sys = ss.lti(A,B,C,D)
        _, sysOut, _ = ss.lsim(sys, U, T)
        ax[i].plot(T, sysOut)

        # discrete:
        sysd = ss.cont2discrete((A,B,C,D), Ts)
        b, a = ss.ss2tf(sysd[0], sysd[1], sysd[2], sysd[3])
        b = b.squeeze() # bug?
        sysOut = list()
        out, state = ss.lfilter(b, a, [0], zi=ss.lfiltic(b,a,[]))
        for j in U:
            out, state = ss.lfilter(b,a, [j], zi=state)
            sysOut.append(float(out))

        ax[i].plot(T, sysOut)
        ax[i].grid(True)

    plt.show(block=False)
