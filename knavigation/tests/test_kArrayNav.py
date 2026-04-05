#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

from knavigation import kArrayNav, kArray
import numpy as np
import math
from numpy import pi, dot
from math  import sqrt


class TestClass_kArrayNav:

    def test_to_rad_to_deg(self):
        print("==== to_rad, to_deg ====")
        euler_deg   = kArrayNav( [-30, 10, 170] )
        euler_rad   = euler_deg.to_rad()
        euler_deg_t = euler_rad.to_deg()

        for i,j in zip(euler_deg, euler_deg_t):
            assert abs(i-j) < 1e-10

    def test_euler_Q_euler(self):
        print("==== euler - Q - euler ====")
        for i in range(20):
            euler_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            euler   = kArrayNav( euler_rad )
            q4      = euler.euler2Q()
            euler_t = q4.Q2euler()
            for j,k in zip(euler, euler_t):
                assert abs(j-k) < 1e-10

    def test_Q2C_C2Q(self):
        print("==== Q2C / C2Q ====")
        for i in range(20):
            euler_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            euler   = kArrayNav( euler_rad )
            q4    = euler.euler2Q()
            C     = q4.Q2C()
            euler_t = C.C2euler()
            q4_t  = C.C2Q()

            for j,k in zip(euler, euler_t):
                assert abs(j-k) < 1e-10

            for j,k in zip(q4, q4_t):
                assert abs(j-k) < 1e-10
    
    def test_q1_x_q2(self):
        print("==== q1_x_q2 ====")
        for i in range(20):
            euler_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            Ca2b = kArrayNav( euler_rad ).euler2C()
            qa2b = kArrayNav( euler_rad).euler2Q()

            euler_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            Cb2c = kArrayNav( euler_rad ).euler2C()
            qb2c = kArrayNav( euler_rad ).euler2Q()

            Ca2c = Cb2c * Ca2b
            qa2c = qb2c.q1_x_q2(qa2b)

            euler   = Ca2c.C2euler()
            euler_t = qa2c.Q2euler()

            for j,k in zip(euler, euler_t):
                assert abs(j-k) < 1e-10

    def test_Re2n(self):
        print("==== Re2n ====")
        Re2n = kArrayNav().Re2n(0,0)
        assert Re2n == kArrayNav( [[0,0,1], [0,1,0], [-1,0,0]] )

    def test_llh_xyz_e(self):
        print("==== llh / xyz_e ====")
        for i in range(20):
            angles_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            llh   = kArrayNav( angles_rad )
            xyz   = llh.ecef_llh2xyz()
            llh_t = xyz.ecef_xyz2llh()
        
            for i,j in zip(llh, llh_t):
                #print("{:f} == {:f} ?".format(i,j))
                assert abs(i-j) < 1e-3

    def test_dynamics(self):
        #----------------------#
        # some dynamic tests:
        #----------------------#
        print("==== dqdt ====")

        # sometimes the developers change the interface with the integrators
        # we shall test it as well:
        import scipy.integrate  as Int

        #  I: inertial frame
        #  b: body frame
        qI2b = kArrayNav( [0,0,0] ).euler2Q()

        # angular rotation between I and b:
        # \omega_{Ib}^I
        w_ib_i = kArray( [2.*pi/180., 0, 0] )

        def eqdiff(q,t,w_ib_i):
            """
            we will use scipy to call this function, and therefore q shall be a list().
            """
            qi2b = kArrayNav(q)
            RI2b = qi2b.Q2C()
            dqdt = qi2b.dqdt( RI2b * kArray( w_ib_i, hvector=False ))
            return dqdt.to_list()

        # a vector described at I:
        F_i = kArray( [0,0,1], hvector=False )
        print("F_i = ")
        print(F_i)

        for t in [1,5,20,89]:
            print()
            # after t seconds, the quaternions should be:
            #print("qI2b.to_list() =")
            #print(qI2b.to_list())
            y = Int.odeint(eqdiff, qI2b.to_list(), [0,t], (w_ib_i,))[1,:]
            # with these euler angles:
            euler = ((180./pi) * kArrayNav(y).Q2euler()).to_list()
            print('euler = {:s}'.format(str(euler)))

            # and described at b:
            F_b = (kArrayNav(y).Q2C() * F_i).to_list()
            euler_expected = kArrayNav((t * w_ib_i)).to_deg() 

            print(f'F_b(phi = {euler[0]:1.03f}) = [ {F_b[0]:1.03f} {F_b[1]:1.03f} {F_b[2]:1.03f} ]')
            print(f'euler_expected = {euler_expected:2.2f}')

            assert (euler_expected - euler).norm() < 1e-3

    def test_gravity(self):
        print("==== gravity ====")

        tmp = kArrayNav([0])
        for lat in kArrayNav( [0, 45, 80] ).to_rad():
            for h in np.linspace(0, 3000, 20):
                tmp = kArrayNav.gravity_n(lat,h)
                assert 9.6 < tmp[2] < 10.0

    #----------------------#
    # coherence tests
    # (sense among frames)
    #----------------------#
    def test_coherence_rotation_body_frame(self):
        # Tests #1

        # this is how the body-frame rotates from the n-frame:
        list_euler = [
                kArrayNav( [0,0,45] ).to_rad(),
                kArrayNav( [0,45,0] ).to_rad(),
                kArrayNav( [45,0,0] ).to_rad(),
        ]

        # with the euler angles above, the transformation of vector [1,0,0]_b will lead to these results:
        sq2 = sqrt(2.)/2.
        list_results_at_n = [
                kArrayNav( [sq2, sq2, 0], hvector=False ),
                kArrayNav( [sq2, 0, -sq2], hvector=False ),
                kArrayNav( [1,0,0], hvector=False ),
        ]

        for euler,vn in zip(list_euler, list_results_at_n):
            # from a fixed reference frame 'n' to 'b'
            Cn2b = euler.euler2C()

            # transforming a vector from 'b' to 'n':
            rb = kArrayNav( [1,0,0], hvector=False )
            rn = Cn2b.T * rb
            assert rn == vn

    def test_coherence_rotation_changes_psi(self):
        # Tests #2
        # An angular velocity [0,0,1] shall increase \psi
        # The angle \psi refers to how 'n' sees 'b'.
        Cn2b   = kArrayNav( [0,0,0] ).to_rad().euler2C()
        Cb2n   = Cn2b.T # to use w_nb_b

        w_nb_b = kArrayNav( [0,0,1], hvector=False )
        Cb2n_p = Cb2n * w_nb_b.to_skew()   # derivative of a tranformation matrix
        Cb2n_step = Cb2n + (0.01 * Cb2n_p) # small integration step
        euler  = Cb2n_step.T.C2euler()     # transposing to obtain the euler again
        assert euler[0][2] > 0
        
    def test_coherence_rotation_changes_all_euler(self):
        # Tests #3
        # An angular velocity [1,1,1] shall increase all euler angles.
        # The euler angles refer to how 'n' sees 'b'.
        Cn2b   = kArrayNav( [0,0,0] ).to_rad().euler2C()
        Cb2n   = Cn2b.T # to use w_nb_b

        w_nb_b = kArrayNav( [1,1,1], hvector=False )
        Cb2n_p = Cb2n * w_nb_b.to_skew()   # derivative of a tranformation matrix
        Cb2n_step = Cb2n + (0.01 * Cb2n_p) # small integration step
        euler  = Cb2n_step.T.C2euler()     # transposing to obtain the euler again
        for i in euler[0]:
            assert i > 0

    def test_coherence_rotation_changes_all_euler_with_quaternions(self):
        # Tests #4
        # Similar to above, but with quaternions in parallel.
        # Note: w_nb is described at b.
        w_nb_b = kArrayNav( [5,5,5], hvector=False ).to_rad() # [rad/s]

        # initial euler angles, from n to b:
        euler_n2b = kArrayNav( [10,10,10] ).to_rad()

        # [k]
        Rb2n = euler_n2b.euler2C().T
        qn2b = euler_n2b.euler2Q()

        # d../dt:
        Rb2n_p = Rb2n * w_nb_b.to_skew() # to use w_nb_b, here the derivative needs Rb2n
        qn2b_p = qn2b.dqdt(w_nb_b)       # to use w_nb_b, here the derivative needs qn2b

        # one small integration step [k+1] (first order):
        dt = 0.01 # [s]
        Rb2n += dt * Rb2n_p
        qn2b += dt * qn2b_p

        # euler:
        euler_k1_R_b2n = Rb2n.T.C2euler().to_deg() # transposing to get Rb2n
        euler_k1_q_n2b = qn2b.Q2euler().to_deg()

        # tests:
        for i in range(3):
            assert euler_k1_R_b2n[0][i] > euler_n2b[0][i]
            assert euler_k1_q_n2b[0][i] > euler_n2b[0][i]
            assert abs( euler_k1_R_b2n[0][i] - euler_k1_q_n2b[0][i] ) < 1e-4

    def test_coherence_inv_quaternion(self):
        # Tests #5
        # Here we want to check whether inverting a quaternion has the same effect as
        # inverting a transformation matrix.
        for i in range(20):
            phi   = 20*np.random.randn()
            theta = 20*np.random.randn()
            psi   = 20*np.random.randn()

            euler = kArrayNav( [phi, theta, psi] ).to_rad()
            R = euler.euler2C()

            q = euler.euler2Q()
            q_inv = q.Qinv()

            euler_from_q_inv = q_inv.Q2euler().to_deg()
            euler_from_R_inv = R.T.C2euler().to_deg()
            assert euler_from_q_inv == euler_from_R_inv

            if False:
                print("euler angles at [k]")
                print(euler.to_deg())
                print("euler angles at [k+1] (from q)")
                print(euler_from_q_inv)
                print("euler angles at [k+1] (from R)")
                print(euler_from_R_inv)
                print()

    def test_coherence_integrating_euler(self):
        # Tests #6
        # To compare the effect of w_nb_b on the euler angles using two methods:
        #  1) dEulerDt()
        #  2) Rb2n_p = Rb2n * w_nb_b
        import scipy.integrate  as Int

        phi   = 20*np.random.randn()
        theta = 20*np.random.randn()
        psi   = 20*np.random.randn()

        def test6_eqdiff_dEulerDt(y, t, *args):
            w = args[0](t)
            dEuler = kArrayNav(y).dEulerDt(w).to_list()
            return dEuler

        def test6_eqdiff_dRdt(y, t, *args):
            w = args[0](t)
            R = kArrayNav( np.asarray(y).reshape(3,3) )
            Rp = R * kArrayNav(w).to_skew()
            return Rp.to_list()

        def tests_angular_speed(t):

            if t < 1:
                w   = kArrayNav( [10,10,10] ).to_rad().to_list()
            elif t < 2:
                w   = kArrayNav( [-10,10,-10] ).to_rad().to_list()
            else:
                w   = kArrayNav( [0,-10,0] ).to_rad().to_list()
            return w

        T   = np.linspace(0, 5, 5*100)

        ret = Int.odeint( test6_eqdiff_dEulerDt, kArrayNav([phi, theta, psi]).to_rad().to_list(), T, args=(tests_angular_speed,) )
        ret = [kArrayNav(i).to_deg().to_list() for i in ret]

        Ret = Int.odeint( test6_eqdiff_dRdt, kArrayNav([phi, theta, psi]).to_rad().euler2C().T.to_list(), T, args=(tests_angular_speed,) )
        Ret = [(kArrayNav( np.asarray(i).reshape(3,3) ).T.C2euler().to_deg() + 0.1).to_list() for i in Ret]
        # this 0.1 above is just to separate the curves on the next graph

        if False:
            plt.figure(2).clf()
            fig,ax = plt.subplots(1,1,num=2)
            ax.plot(T, ret)
            ax.plot(T, Ret)
            ax.grid(True)

        for i,j in zip(ret, Ret):
            for m,n in zip(i,j):
                assert abs(m-n) < 0.11

        ###########################
        #plt.show(block=False)
        ###########################

    def test_polar_projection(self):
        print("==== polar projection ====")
        phi = 10
        tta = 20
        C   = lambda psi: kArrayNav( [phi, tta, psi] ).to_rad().euler2C()

        psi_0 = 30
        dpsi  = 1
        dC    = (C(psi_0 + dpsi) - C(psi_0)) / dpsi

        psi_1 = 60
        C_with_psi_1 = C(psi_0) + ((psi_1-psi_0) * dC) # using euler to integrate up to psi_1
        aux = C_with_psi_1 * C_with_psi_1.T
        print("before orthogonalization, C x C^T = \n", aux)
        print("euler angles = ", C_with_psi_1.C2euler().to_deg())

        # this test is valid if the main diagonal is of numbers far from 1.0:
        for i in range(3):
            for j in range(3):
                if i == j:
                    assert abs(aux[i][j] - 1.0) > 0.03

        C90 = C_with_psi_1.polar_projection() # orthogonalizing the transformation matrix
        aux = C90 * C90.T
        print("after  orthogonalization, C x C^T = \n", aux)
        print("euler angles = ", C90.C2euler().to_deg())

        # we are looking for an identity matrix:
        for i in range(3):
            for j in range(3):
                if i == j:
                    assert abs(aux[i][j] - 1.0) < 1e-10
                else:
                    assert abs(aux[i][j]) < 1e-10

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
