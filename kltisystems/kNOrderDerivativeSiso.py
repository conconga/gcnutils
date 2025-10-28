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
class kSystemBA:
    def __init__(self, num, den):
        self.B = num
        self.A = den

        # initialize the filter:
        _, self.state = ss.lfilter(self.B, self.A, [0], zi=ss.lfiltic(self.B, self.A, []))

    def update(self, sample):
        out, self.state = ss.lfilter(self.B, self.A, [sample], zi=self.state)
        return(out)

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

        self.order = order
        self.pole  = pole
        self.Ts    = Ts

        self._calculate_stable_transfer_function()
        self.A, self.B, self.D = self._calculate_state_space_matrices()
        self._setup_systems_for_each_derivative()

    def _calculate_stable_transfer_function(self):
        # polynomial:
        den = np.poly1d( [1, -pole] )
        for i in range(self.order):
            # den = [1 ..  a1 a0]
            den *= np.poly1d( [1, -pole] )
        self.den    = list(den) # indexing 'den' is non intuitive
        self.revden = self.den[::-1] # reversing

    def _calculate_state_space_matrices(self):
        A     = np.hstack( (np.zeros((self.order, 1)), np.eye(self.order)) )
        A     = np.vstack( (A, [-i for i in self.revden[:self.order+1]]) )
        B     = np.zeros(self.order+1)
        B[-1] = self.den[-1]
        B     = B.reshape((-1,1))
        D     = 0;

        return A, B, D

    def _setup_systems_for_each_derivative(self):
        system = list()
        for i in range(self.order+1):
            # when i==0, the system will be a LP filter to smooth the input
            # signal and enable the calculation of the derivatives.
            C = np.zeros(order+1)
            C[i] = 1.0
            
            # setup of discrete systems (filters):
            sysd = ss.cont2discrete(( self.A, self.B, C, self.D ), self.Ts)
            b, a = ss.ss2tf(sysd[0], sysd[1], sysd[2], sysd[3])
            b = b.squeeze() # bug?
            system.append( kSystemBA( b,a ) )

        self.system = system

    def update(self, sample):
        """
        return:
            out[0] = smoothed input
            out[1] = first order derivative
            out[2] = second order derivative
            ...
            out[N] = Nth-order derivative
        """
        out = list()
        for sys in self.system:
            out.append( sys.update(sample) )

        return out

#################################
## ##WWww=--  main:  --=wwWW## ##
#################################
class kNOrderDerivativeSisoTests:
    def do_tests(self):

        import matplotlib.pylab as plt

        pole  = -100
        order = 4
        Fs    = 400 # [Hz]
        Ts    = 1./Fs
        obj   = kNOrderDerivativeSiso( 4, -100, Ts )

        # time and input signal:
        U = list()
        for i in range(4):
            v = np.random.randn()
            U.append( [v for i in range(Fs)] )
        U = np.hstack(U)
        T = np.asarray(range(len(U))) * Ts

        plt.figure(1).clf()
        fig,ax = plt.subplots(order+1, 1, num=1, sharex=True)
        ax[0].plot(T, U)

        # all derivatives:
        sysOut = list()
        for i in U:
            sysOut.append( obj.update(i) )
        sysOut = np.asarray(sysOut)

        # for the pictures:
        for i in range(order+1):
            ax[i].plot(T, sysOut[:,i])
            ax[i].grid(True)

        plt.show(block=False)
