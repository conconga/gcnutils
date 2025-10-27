#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung: An object to simplify steps to integrate a SoSODE object.
Autor: Luciano Auguto Kruk
Erstellt am: 26.10.2025
Version: 1.0.0
Lizenz: Please keep this header with the file.
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##

import numpy            as np
import scipy.integrate  as Int
from kSosode import *

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSosodeUtils:
    """
    Helpful methods to be adopted when modeling a system with SoSODE.
    """

    def __init__(self, *args, **kargs):
        super().__init__(*args, **kargs)

    def _try_to_get(self, kargs, txt):
        try:
            ret = kargs[txt]
        except:
            ret = False

        return ret

    def pick_from_state(self, *args):
        """
        Returns the current value for a particular state.

        For example: 
            `self.pick_from_state('vN')`          returns the most updated value of 'vN'.
            `self.pick_from_state('vN', 'vN_p')`  returns a list with the updated respective values.

        Condition:
            The list `self.order_states` shall be available with the ordered
            names of the variables in the state vector.  Moreover, the state
            vector shall obey this order.
        """

        if len(args) == 0:
            raise(NameError("the method needs at least one argument"))

        elif len(args) == 1:
            idx = self.order_states.index(args[0])
            ret = self.state[idx]

        else:
            ret = list()
            for txt in args:
                idx = self.order_states.index(txt)
                ret.append(self.state[idx])

        return ret

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSosodeIntegrator:
    def __init__(self, *args, **kargs):
        super().__init__(*args, **kargs)

        self.curr_time = -1.0 # it shall start negative

    def update(self):
        """
        Step update at the integration.
            self.dt         :  integration step
            self.curr_time  :  time instant of the last available integration step
                            :  the curr_time shall be initialized negative before any integration
            self.state0     :  provided state vector for t=0[s]
        """
        if self.curr_time < 0:
            # only for t=0
            self.curr_time = 0.0
            self.state = self.state0
            return self.state

        # target-time, one step:
        t = self.curr_time + self.dt

        # integrate one step:
        R = Int.odeint( self.sys, self.state, [self.curr_time, t], args=() )
        self.curr_time = t
        self.state     = R[1]

        return self.state

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kExample_RC_discharge_system:
    """
    Modeling of the following differential equation:

        d v(t) / dt = -(1/RC).v(t)

    """
    def __init__(self, *args, **kargs):
        super().__init__(*args, **kargs)

    def __init__(self, **kargs):
        R   = 100 if not self._try_to_get(kargs, "R")    else self._try_to_get(kargs, "R")
        C   = 100e-9 if not self._try_to_get(kargs, "C") else self._try_to_get(kargs, "C")
        V   = 5.0 if not self._try_to_get(kargs, "V")    else self._try_to_get(kargs, "V")

        self.R = R
        self.C = C

        self.order_states = ["V"]
        self.state0 =       [ V ]

        # system registration:
        fn_dVdt = kSosodeFunction(self.sys_dVdt)
        fn_dVdt.set_i_state([ 'V' ]) # receives 'V'
        fn_dVdt.set_o_state([ 'V' ]) # returns the first derivative of 'V'

        self.sys = kSosode( fn_dVdt, reverse=True, order_states=self.order_states )
        self.sys.create_nets()

    def sys_dVdt(self, t, state):
        V    = state[0]
        dVdt = -V/(self.R*self.C)
        return [dVdt]

class kExample_RC_discharge(kSosodeUtils, kSosodeIntegrator, kExample_RC_discharge_system):
    def __init__(self, sample_freq_Hz, **kargs):
        # initialize the System of Systems model:
        super().__init__(**kargs)

        assert sample_freq_Hz is not None
        assert sample_freq_Hz > 0

        self.dt        = 1./sample_freq_Hz

    def get_V(self):
        return self.pick_from_state( 'V' )

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
if __name__ == "__main__":
    sample_freq_Hz    = 100
    t_max_simu        = 2
    T                 = np.arange(0, t_max_simu, 1./sample_freq_Hz)

    ###########################
    ##mm==-- main loop --==mm##

    import matplotlib.pylab as plt
    plt.figure(1).clf()
    fig, ax = plt.subplots(1,1,num=1)

    for rc in [kExample_RC_discharge(100, R=1e6), kExample_RC_discharge(100, V=3, R=100e3, C=4.7e-6)]:
        rc_log = list()
        for t in T:
            # update current state vector:
            rc_log.append( rc.update() )

        ax.plot(T, rc_log)

    ax.grid(True)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("Vc [V]")

    #####################
    plt.tight_layout()
    plt.show(block=False)
    #####################

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
