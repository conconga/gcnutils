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
import numpy           as np
#from   numpy           import dot
#from   numpy           import inf
#from   scipy.integrate import odeint
import scipy.signal     as ss
from scipy.signal       import ss2tf
import matplotlib.pylab as plt

#################################
## ##WWww=--  main:  --=wwWW## ##
#################################
if __name__ == "__main__":
    pole = -100
    order = 4
    Fs = 100 # [Hz]
    Ts = 1./Fs
    U = list()

    for i in range(4):
        v = np.random.randn()
        U.append( [v for i in range(Fs)] )
    U = np.hstack(U)
    T = np.asarray(range(len(U))) * Ts

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
        ax[i].grid(True)

    plt.show(block=False)
