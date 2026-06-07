#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

from kltisystems import k2OrderLTIsysMimoFactory
import numpy   as np
from   numpy   import inf, pi
import pytest

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
def fn_set_trace():
    import pudb
    pudb.set_trace()

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class TestClass_Example:
    block = 0

    def test_example_2order_mimo(self):

        from   scipy.integrate    import odeint
        import matplotlib.pyplot  as plt

        Fs  = 200 # [Hz]
        Ts  = 1./Fs
        T   = np.arange(0,2,Ts)
        U   = (T > 0.5) * 1.0

        #UmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUmUm#
        # test of MIMO systems:

        mimo_c = k2OrderLTIsysMimoFactory(0.7, 2.*np.pi*20, np.asarray([0,1]), -5, 5, -2, [0.8,1.1])
        mimo_d = k2OrderLTIsysMimoFactory(0.7, 2.*np.pi*20, np.asarray([0,1]), -5, 5, -2, [0.8,1.1], Ts=Ts)
        mimo_c_buf = [];
        mimo_d_buf = [];

        for t,u in zip(T,U):

            # integration:
            # continuous:
            mimo_c.update(t, u)

            # discrete:
            mimo_d.update(u)

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

        plt.show(block=self.block)



#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>                                                      >>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

class TestClass_2order_mimo:


    pass


#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
