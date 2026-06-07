#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung: These classes implement the dynamical behaviour of a second order
    LTI system SISO and saturated in input and output rate.
Autor: Luciano Auguto Kruk
Erstellt am: 27.09.2025
Version: 1.0.0
Lizenz: Please keep this header with the file.
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##
import numpy           as np
from   numpy           import dot
from   scipy.integrate import odeint

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kCommon2OrderLTIsysSiso:
    """
    parameters:                     
        qsi : damping factor        
        wn  : [rad/s] natural freq.
        x0  : initial state       
        min_dxdt : min rate slope
        max_dxdt : max rate slope
        min_x    : min output   
        max_x    : max output  
        Ts       : discrete interval
    """

    def __init__(self, qsi, wn, x0, min_dxdt, max_dxdt, min_x, max_x, Ts=0):
        """
            Scalar SISO second order LTI filter.

            For continuous simulation, use these methods:
                .dstate_dt() to calculate the derivative, and
                .update() to update the current state.

            For discrete simulation (Ts>0), use this method:
                .update() to update the current state for a given reference input value.

            The parameters 'qsi,wn,min_dxdt, max_dxdt, min_x, max_x' shall be scalars.
        """

        # this class shall not be instanciated.
        assert "kCommon2OrderLTIsysSiso" not in str(self.__class__)

        if isinstance(x0, (int, float, np.int64)):
            x0 = np.asarray([x0, 0])

        if (len(x0) == 1):
            self.x1  = x0[0] # output x
            self.x2  = 0.0   # dx/dt
        else:
            self.x1  = x0[0] # output x
            self.x2  = x0[1] # dx/dt

        self.x   = np.asarray(x0)

        self.a = wn**2.0
        self.b = 2.*qsi*wn
        self.u = 0.0

        self.min_dxdt = min_dxdt
        self.max_dxdt = max_dxdt
        self.min_x    = min_x
        self.max_x    = max_x
        self.Ts       = Ts

    def __str__(self):
        return "state = %10.3e, derivative = %10.3e" % (self.x1, self.x2)

    def _saturate(self, _in, _min, _max):
        if (_in < _min):
            _out = _min
        elif (_in > _max):
            _out = _max
        else:
            _out = _in

        return _out

    def get_state(self):
        return (self.x)


#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class k2OrderLTIsysSisoContinuous (kCommon2OrderLTIsysSiso):
    def __init__(self, *args, **kargs):
        super().__init__(*args, **kargs)
        self.curr_t = 0

    def _dstate_dt(self, x, t, u):
        """
            The state vector is

                x = [x1, x2]^T

            and the system is:

                dot(x1) = x2
                dot(x2) = (....)
        """

        # states:
        v1 = x[0]
        v2 = x[1]

        # saturation:
        f1 = self._saturate(u, self.min_x, self.max_x)
        f2 = v1 + (v2*self.b/self.a)

        vp1 = v2

        if ((v2 <= self.min_dxdt) and (f1 <= f2)) or ((v2 >= self.max_dxdt) and (f1 >= f2)):
            vp2 = 0
        else:
            vp2 = (-self.b*v2) + (self.a*(f1-v1))

        return np.asarray([vp1, vp2])

    def update(self, t, u):
        """
        Continuous time update.

        inputs:
        x: state at time t
        u: reference input
        """

        if t > self.curr_t:
            t0 = self.curr_t
            _,y = odeint(self._dstate_dt, self.get_state(), [t0, t], (u,))
            self.curr_t = t

            self.x  = np.asarray(y)
            self.x1 = y[0]
            self.x2 = y[1]

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class k2OrderLTIsysSisoDiscrete (kCommon2OrderLTIsysSiso):
    def __init__(self, *args, **kargs):
        super().__init__(*args, **kargs)

        if self.Ts > 0:
            a        = self.a
            b        = self.b
            self.Ac  = np.asarray([[0,1],[-a,-b]])
            self.Bc  = np.asarray([0,a])

            # bilinear:
            aux      = np.linalg.inv(np.eye(2) - (0.5*self.Ts*self.Ac))
            self.Ad  = dot(np.eye(2) + (0.5*self.Ts*self.Ac), aux)
            self.Bd  = dot(aux, self.Bc) * self.Ts

    def update(self, t, u):
        """
        Discrete time update.
        u: input at time t
        """

        assert(self.Ts > 0)

        u_sat = self._saturate(u, self.min_x, self.max_x)

        w2k1   = dot(self.Ad, self.x) + (self.Bd * u_sat) # vector
        v2k    = self.x[1] # escalar v2[k]
        v2k1   = w2k1[1]   # escalar v2[k+1]

        w2k1_c = dot(self.Ac, self.x) + (self.Bc * u_sat)

        # lower constrained at [k+1]
        if (v2k1 < self.min_dxdt):

            # and trying to sink:
            if (v2k > v2k1):

                # forward:
                ti = (1./self.Ts) * (self.min_dxdt-v2k) / w2k1_c[1]

                if (0.<ti<1.): # descendo, com saturacao no caminho;
                    # forward:
                    Adti = np.eye(2) + (self.Ts*ti*self.Ac)
                    Bdti = self.Bc*self.Ts*ti

                    # from 0 until ti*Ts:
                    w2k1_ = dot(Adti, self.x) + (Bdti*u_sat)

                    # from ti*Ts until Ts:
                    v1k1 = w2k1_[0] + ((1.0 - ti) * self.Ts * w2k1_[1])

                    # update:
                    w2k1 = np.asarray([v1k1, w2k1_[1]])

                elif (ti<=0.): # comecou [k] em saturacao;
                    # integrate with d(v2k)/dt = 0:
                    v1k1 = self.x[0] + (self.Ts * self.min_dxdt)
                    w2k1 = np.asarray([v1k1, self.min_dxdt])

                else: # nao satura ateh [k+1]; impossible!
                    # w2k1 <- w2k1
                    pass

            else: # (v2k1 > v2k) # na regiao proibida de saturacao em [k] e em [k+1];
                # integrate with d(v2k)/dt = 0:
                v1k1 = self.x[0] + (self.Ts * self.min_dxdt)
                w2k1 = np.asarray([v1k1, self.min_dxdt])

        # upper constrained at [k+1]:
        elif (v2k1 > self.max_dxdt):

            # and trying to climb:
            if (v2k1 > v2k):

                # forward:
                ti = (1./self.Ts) * (self.max_dxdt-v2k) / w2k1_c[1]

                if (0.<ti<1.): # subindo com saturacao no caminho:
                    # forward:
                    Adti = np.eye(2) + (self.Ts*ti*self.Ac)
                    Bdti = self.Bc*self.Ts*ti

                    # from 0 until ti*Ts:
                    w2k1_ = dot(Adti, self.x) + (Bdti*u_sat)

                    # from ti*Ts until Ts:
                    v1k1 = w2k1_[0] + ((1.0 - ti) * self.Ts * w2k1_[1])

                    # update:
                    w2k1 = np.asarray([v1k1, w2k1_[1]])

                elif (ti<=0.): # comecou [k] em saturacao;
                    # integrate with d(v2k)/dt = 0:
                    v1k1 = self.x[0] + (self.Ts * v2k)
                    w2k1 = np.asarray([v1k1, v2k])

                else: # nao satura ateh [k+1]; impossible!
                    # w2k1 <- w2k1
                    pass

            else: # (v2k > v2k1) # na regiao proibida de saturacao em [k] e em [k+1];
                # integrate with d(v2k)/dt = 0:
                v1k1 = self.x[0] + (self.Ts * self.max_dxdt)
                w2k1 = np.asarray([v1k1, self.max_dxdt])

        else: # easy, no constraint:
            pass


        self.x = w2k1
        self.x1 = self.x[0]
        self.x2 = self.x[1]

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#

def k2OrderLTIsysSisoFactory(*args, **kargs):

    Ts = kargs.get("Ts", None)
    if (Ts is None):
        return k2OrderLTIsysSisoContinuous(*args, **kargs)
    else:
        if Ts > 0:
            return k2OrderLTIsysSisoDiscrete(*args, **kargs)
        else:
            return k2OrderLTIsysSisoContinuous(*args, **kargs)

#====================================#
