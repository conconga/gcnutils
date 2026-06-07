#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

from kltisystems import k1OrderLTIsysSisoFactory
from kltisystems import k1OrderLTIsysMimoDiscrete
from knavigation import kArray
from math        import exp
import numpy as np
import pytest

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class TestClass_Example:
    block = 0

    def test_example_siso(self):

        import matplotlib.pyplot      as plt

        a    = -10.  # pole for the transfer function:
        tmax = 2.0
        Ts   = 5e-3 # sample rate
        T    = [i*Ts for i in range(int(tmax/Ts)+1)] # time vector
        f1c  = k1OrderLTIsysSisoFactory(a, 0)
        f1d  = k1OrderLTIsysSisoFactory(a, Ts, 1)

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
        plt.show(block=self.block)
        #.............................................#

    def test_example_mimo(self):

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
        plt.show(block=self.block)
        #.............................................#


#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
def fn_solution(y0, x, pole, t):
    return x + ((y0-x) * exp(pole * t))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class TestClass_1OrderContinuos:

    def test_decay(self):
        pole = -2.0
        y0   =  3.0

        f1c = k1OrderLTIsysSisoFactory(pole, y0)

        for t in [0, 0.5, 1, 2]:
            assert abs(f1c.update(t, 1.0) - fn_solution(y0, 1.0, pole, t)) < 1e-7

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class TestClass_1OrderDiscrete:

    def test_decay(self):
        pole = -2.0
        y0   = 3.0
        Ts   = 1./1e6 # [s]

        f1d = k1OrderLTIsysSisoFactory(pole, Ts, y0)

        for t in np.arange(Ts, 2.0, Ts):
            y = f1d.update(1.0)
            #print(f"y = {y};    sol = {fn_solution(y0, 1.0, pole, t)}")
            assert abs(y - fn_solution(y0, 1.0, pole, t)) < 1e-6

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class TestClass_1OrderDiscreteMimo:

    def test_decay(self):
        pole = -2.0
        y0   = [3.0, 5.0]
        Ts   = 1./1e6 # [s]

        f1d = k1OrderLTIsysMimoDiscrete(pole, Ts, y0)

        for t in np.arange(Ts, 0.5, Ts):
            y = f1d.update([0,0])

            #print(f"y = {y};    sol = {fn_solution(y0, 1.0, pole, t)}")
            for i in range(2):
                assert abs(y[i] - fn_solution(y0[i], 0.0, pole, t)) < 1e-6

    def test_input_nparray(self):
        pole = -2.0
        y0   = [3.0, 5.0]
        Ts   = 1./1e3 # [s]

        f1d = k1OrderLTIsysMimoDiscrete(pole, Ts, y0)

        for t in np.arange(Ts, 1.0, Ts):
            y = f1d.update(np.asarray([0,0]))

            #print(f"y = {y};    sol = {fn_solution(y0, 1.0, pole, t)}")
            for i in range(2):
                assert abs(y[i] - fn_solution(y0[i], 0.0, pole, t)) < 1e-6

    def test_input_karray(self):
        pole = -2.0
        y0   = [3.0, 5.0]
        Ts   = 1./1e3 # [s]

        f1d = k1OrderLTIsysMimoDiscrete(pole, Ts, y0)

        for t in np.arange(Ts, 1.0, Ts):
            y = f1d.update(kArray([0,0], hvector=False))

            #print(f"y = {y};    sol = {fn_solution(y0, 1.0, pole, t)}")
            for i in range(2):
                assert abs(y[i] - fn_solution(y0[i], 0.0, pole, t)) < 1e-6

    def test_input_multiple_poles(self):
        pole = [-2.0, -1]
        y0   = [3.0, 5.0]
        Ts   = 1./1e3 # [s]

        f1d = k1OrderLTIsysMimoDiscrete(pole, Ts, y0)

        for t in np.arange(Ts, 1.0, Ts):
            y = f1d.update(kArray([0,0], hvector=False))

            #print(f"y = {y};    sol = {fn_solution(y0, 1.0, pole, t)}")
            for i in range(2):
                assert abs(y[i] - fn_solution(y0[i], 0.0, pole[i], t)) < 1e-6

    def test_input_multiple_poles_one_instable(self):
        pole = [-2.0, 2]
        y0   = [3.0, 5.0]
        Ts   = 1./1e3 # [s]

        with pytest.raises(Exception) as e:
            f1d = k1OrderLTIsysMimoDiscrete(pole, Ts, y0)

    def test_single_state_vector(self):
        pole = [-2.0, -1]
        y0   = [3.1, 5.1]
        Ts   = 1./1e3 # [s]

        f1d = k1OrderLTIsysMimoDiscrete(pole, Ts, y0)

        for i in range(2):
            assert int(f1d.y[i]) == int(y0[i])


#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
