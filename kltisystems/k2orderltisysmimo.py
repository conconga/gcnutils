#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei:
Beschreibung: These classes implement the dynamical behaviour of a second order
    LTI system MIMO and saturated in input and output rate.
Autor: Luciano Auguto Kruk
Erstellt am: 27.09.2025
Version: 1.0.0
Lizenz: Please keep this header with the file.
GitHub:
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##
import numpy   as np
from   numpy   import inf
from   k2orderltisyssiso import k2OrderLTIsysSiso

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kCommon2OrderLTIsysMimo:

    def __init__(self, qsi, wn, x0, min_dxdt, max_dxdt, min_x, max_x, Ts=0):
        """
        \parameters:
            qsi : damping factor
            wn  : [rad/s] natural freq.
            x0  : initial state
            min_dxdt : min rate slope
            max_dxdt : max rate slope
            min_x    : min output
            max_x    : max output
            Ts       : discrete interval
        """

        # this class shall not be instanciated.
        assert "kCommon2OrderLTIsysMimo" not in str(self.__class__)

        self.n        = len(x0)
        self.x0       = x0
        self.Ts       = Ts
        self.qsi      = self._fn_fill_config_list(qsi)
        self.wn       = self._fn_fill_config_list(wn)
        self.min_x    = self._fn_fill_config_list(min_x)
        self.max_x    = self._fn_fill_config_list(max_x)
        self.min_dxdt = self._fn_fill_config_list(min_dxdt)
        self.max_dxdt = self._fn_fill_config_list(max_dxdt)

        # instances of k2OrderLTIsysSiso:
        self.siso = [ k2OrderLTIsysSiso(
            self.qsi[i],
            self.wn[i],
            self.x0[i],
            self.min_dxdt[i],
            self.max_dxdt[i],
            self.min_x[i],
            self.max_x[i],
            Ts=Ts) for i in range(self.n)
        ]

    def _fn_fill_config_list(self, In):
        """
            Converts scalar, array or list into list.
        """

        if (isinstance(In, int) or isinstance(In, float)):
            Out = [In for i in range(self.n)]

        elif (isinstance(In, np.ndarray)):
            Out = In.tolist()

        elif (isinstance(In, list)):
            Out = In

        else:
            raise(NameError("hummm... are you sure you know what you are doing? I'm not!"))

        return Out

    def get_state(self,*args):
        """
            Returns state of MIMO object or individual SISO instances.

            .get_state() returns state from MIMO object.
            .get_state(i) returns state from SISO instance 'i'.
        """

        if (len(args) == 0):
            x = []
            for i in range(self.n):
                x += self.siso[i].get_state().tolist()
        elif (len(args) == 1):
            x = self.siso[args[0]].get_state()

        return np.asarray(x) # list or np.array

    def interleave(self, x):
        """
            This method changes
                x = [x1, x2, ... , x1p, x2p, ...]
            to
                x = [x1, x1p, x2, x2p, ... ]
        """

        if (isinstance(x, np.ndarray)):
            In = x.squeeze().tolist()

        # interleave:
        y = []
        for i in range(self.n):
            y += [In[i], In[self.n+i]]

        if (isinstance(x, np.ndarray)):
            y = np.asarray(y)

        return y

    def deinterleave(self, x):
        """
            This method changes
                x = [x1, x1p, x2, x2p, ... ]
            to
                x = [x1, x2, ... , x1p, x2p, ...]
        """

        if (isinstance(x, np.ndarray)):
            In = x.squeeze().tolist()

        # deinterleave:
        y = [ In[2*i] for i in range(self.n) ] + [ In[(2*i)+1] for i in range(self.n) ]

        if (isinstance(x, np.ndarray)):
            y = np.asarray(y)

        return y

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class k2OrderLTIsysMimoContinuos (kCommon2OrderLTIsysMimo):
    def __init__(self, *args, **kargs):
        super().__init__(*args, **kargs)
        self.curr_t = 0

    def _c_update(self, t, u):
        """
        Continuous time update.
        x: state at time t
        """

        # protecting 'u':
        u = self._fn_fill_config_list(u)

        # all updates:
        for i in range(self.n):
            self.siso[i].update(t, u[i])

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class k2OrderLTIsysMimoDiscrete (kCommon2OrderLTIsysMimo):
    def __init__(self, *args, **kargs):
        super().__init__(*args, **kargs)

    def _d_update(self, t, u):
        """
        Discrete time update.
        u: input at time t
        """

        assert(self.Ts > 0)

        # protecting 'u':
        u = self._fn_fill_config_list(u)

        for i in range(self.n):
            self.siso[i].update(t, u[i])

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class k2OrderLTIsysMimo (k2OrderLTIsysMimoContinuos, k2OrderLTIsysMimoDiscrete):
    def __init__(self, *args, **kargs):
        """
            MIMO second order LTI filter.

            For continuous simulation, use these methods:
                ._c_update() to update the current state.

            For discrete simulation (Ts>0), use this method:
                ._d_update() to update the current state for a given reference input value.

            The parameters 'qsi,wn,min_dxdt, max_dxdt, min_x, max_x' shall be
            either scalars, vectors or lists. The parameter 'Ts' shall be a scalar, if any.

            The dimension of 'x0' defines the number of internal SISO systems.

            The MIMO state vector represents each two values one single SISO state vector:

            state =   [   x1    ]
                      [ dot{x1} ]
                      [   x2    ]
                      [ dot{x2} ]
                      [   ...   ]
                      [   xn    ]
                      [ dot{xn} ]

            The methods interleave() and deinterleave() change the representation of state.
        """
        super().__init__(*args, **kargs)

    def update(self, *args):

        if self.Ts > 0:
            # discrete:
            self._d_update(*args)
        else:
            # continuous:
            self._c_update(*args)

#################################
## ##WWww=--  main:  --=wwWW## ##
#################################
if (__name__ == "__main__"):

    from   scipy.integrate    import odeint
    import matplotlib.pyplot  as plt

    Fs  = 200 # [Hz]
    Ts  = 1./Fs
    T   = np.arange(0,2,Ts)
    U   = (T > 0.5) * 1.0

    #UmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUm#
    # test of MIMO systems:

    mimo_c = k2OrderLTIsysMimo(0.7, 2.*np.pi*20, np.asarray([0,1]), -5, 5, -2, [0.8,1.1])
    mimo_d = k2OrderLTIsysMimo(0.7, 2.*np.pi*20, np.asarray([0,1]), -5, 5, -2, [0.8,1.1], Ts=Ts)
    mimo_c_buf = [];
    mimo_d_buf = [];

    for t,u in zip(T,U):

        # integration:
        # continuous:
        mimo_c.update(t, u)

        # discrete:
        mimo_d.update(t, u)

        # buffer:
        mimo_c_buf.append([t,u] + mimo_c.get_state().tolist())
        mimo_d_buf.append([t,u] + mimo_d.get_state().tolist())

    mimo_c_buf = np.asarray(mimo_c_buf)
    mimo_d_buf = np.asarray(mimo_d_buf)

    #  - #  - #  - #  - #  - #  - #  - #  - # - #
    # figures:
    #  - #  - #  - #  - #  - #  - #  - #  - # - #

    #----- new figure -----#
    plt.figure(3), plt.clf()
    plt.plot(mimo_c_buf[:,0], mimo_c_buf[:,1:])
    plt.grid(True)
    plt.legend(('signal', 'c1', 'dot{c1}', 'c2', 'dot{c2}'))

    #----- new figure -----#
    plt.figure(4), plt.clf()
    plt.plot(mimo_d_buf[:,0], mimo_d_buf[:,1:])
    plt.grid(True)
    plt.legend(('signal', 'd1', 'dot{d1}', 'd2', 'dot{d2}'))

    plt.show(block=False)

#====================================#
