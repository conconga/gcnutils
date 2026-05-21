#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

from knavigation import kArrayNav, kQuatNav
import numpy as np
import math
import pytest
from unittest.mock import patch
from numpy         import pi, dot
from math          import sqrt


class TestClass_kQuatNav:

    def get_random_euler_rad(self):
        # vector [3] with random euler angles [-pi/2, pi/2]
        #return np.asarray([10,20,30]) * pi/180
        return (((np.random.rand(1,3)*2.0)-1)*90) * pi/180

    def test_quat_conjugate(self):
        q  = kArrayNav( [1,2,3,4], hvector=False )
        qj = q.q_conj()

        for i in range(4):
            if i == 0:
                assert abs(q[i] - qj[i]) < 1e-10
            else:
                assert abs(q[i] + qj[i]) < 1e-10

    def test_quat_norm(self):
        for i in range(20):
            euler   = kArrayNav( self.get_random_euler_rad() )
            q4      = euler.euler2Q() # <= this is normalized
            scale   = np.random.randint(2,20)
            q_scale = q4 * scale
            q_norm  = q4.q_norm()
            for j,k in zip(q4, q_norm):
                assert abs(j-k) < 1e-10

    def test_q_x_3d(self):
        vector = kArrayNav( [1,2,-3], hvector=False )
        for _ in range(30):
            euler = kArrayNav( self.get_random_euler_rad() )

            v1 = euler.euler2C() * vector
            v2 = euler.euler2Q().q_x_3d(vector)

            norm2_v1 = sum([i**2 for i in v1])
            norm2_v2 = sum([i**2 for i in v2])

            assert abs(norm2_v1 - norm2_v2) < 1e-10
            assert v1 == v2
    
    def test_q_x_q_against_vector3D(self):
        vector = kArrayNav( [1,-2,3], hvector=False )
        for _ in range(30):
            euler_a2b = kArrayNav( self.get_random_euler_rad() )
            euler_b2c = kArrayNav( self.get_random_euler_rad() )

            q_a2b = euler_a2b.euler2Q()
            q_b2c = euler_b2c.euler2Q()

            # ground truth:
            if np.random.rand() > 0.5:
                vc_1 = euler_b2c.euler2C() * euler_a2b.euler2C() * vector
            else:
                vc_1 = q_b2c.Q2C() * q_a2b.Q2C() * vector

            # with two consecutive transformations:
            vc_2 = q_b2c.q_x_3d( q_a2b.q_x_3d(vector) )

            # now cascading the quaternions:
            vc_3 = q_a2b.q_x_q(q_b2c).q_x_3d(vector)

            # verifying the norm of the vectors:
            norm2_v1 = sum([i**2 for i in vc_1])
            norm2_v2 = sum([i**2 for i in vc_2])
            norm2_v3 = sum([i**2 for i in vc_3])

            assert abs(norm2_v1 - norm2_v2) < 1e-10
            assert abs(norm2_v1 - norm2_v3) < 1e-10

            # comparing the components in the vectors:
            for i,j in zip(vc_1, vc_2):
                assert abs(i-j) < 1e-10

            for i,j in zip(vc_1, vc_3):
                assert abs(i-j) < 1e-10

    def test_q_x_q_tests_against_C(self):
        for _ in range(50):
            euler_a2b = kArrayNav( self.get_random_euler_rad() )
            Ca2b      = kArrayNav( euler_a2b ).euler2C()
            qa2b      = kArrayNav( euler_a2b ).euler2Q()

            euler_b2c = kArrayNav( self.get_random_euler_rad() )
            Cb2c      = kArrayNav( euler_b2c ).euler2C()
            qb2c      = kArrayNav( euler_b2c ).euler2Q()

            Ca2c = Cb2c * Ca2b       # <= truth
            qa2c = qa2b.q_x_q(qb2c)  # <= test

            C    = qa2c.Q2C()

            for j,k in zip(Ca2c, C):
                assert abs(j-k) < 1e-10

    def test_q_x_q_tests_against_euler(self):
        for _ in range(50):
            euler_a2b = kArrayNav( self.get_random_euler_rad() )
            Ca2b = kArrayNav( euler_a2b ).euler2C()
            qa2b = kArrayNav( euler_a2b ).euler2Q()

            euler_b2c = kArrayNav( self.get_random_euler_rad() )
            Cb2c = kArrayNav( euler_b2c ).euler2C()
            qb2c = kArrayNav( euler_b2c ).euler2Q()

            Ca2c = Cb2c * Ca2b
            qa2c = qa2b.q_x_q(qb2c)

            euler   = Ca2c.C2euler()
            euler_t = qa2c.Q2euler()

            for j,k in zip(euler, euler_t):
                assert abs(j-k) < 1e-10


    def test_q_x_q_calculated_in_two_ways(self):
        for _ in range(50):
            euler_rad = kArrayNav( self.get_random_euler_rad() )
            qa2b = kArrayNav( euler_rad).euler2Q()

            euler_rad = kArrayNav( self.get_random_euler_rad() )
            qb2c = kArrayNav( euler_rad ).euler2Q()

            qa2c = qa2b.q_x_q(qb2c)
            qnew = qa2b._q_x_q(qb2c)

            for i,j in zip(qa2c, qnew):
                assert abs(i-j) < 1e-10

            e1 = qa2c.Q2euler().to_deg()
            e2 = qnew.Q2euler().to_deg()

            for i,j in zip(e1,e2):
                assert abs(i-j) < 1e-10

    def test_q_x_q1(self):
        for _ in range(100):
            euler_rad = kArrayNav( self.get_random_euler_rad() )
            q  = euler_rad.euler2Q()

            if np.random.rand() > 0.5:
                qj = q.q_conj() # it works only for quaternions of transformation
            else:
                qj = q.q_inv()

            prod1 = q.q_x_q(qj)
            prod2 = qj.q_x_q(q)
            ret_expected = kArrayNav([1,0,0,0], hvector=False)

            #print(f"expected == {prod1}   ==>  {ret_expected == prod1}")
            #print(f"expected == {prod2}   ==>  {ret_expected == prod2}")

            for i in range(4):
                if i == 0:
                    assert abs(prod1.squeeze()[i] - 1.0) < 1e-10
                    assert abs(prod2.squeeze()[i] - 1.0) < 1e-10
                else:
                    assert abs(prod1.squeeze()[i]) < 1e-10
                    assert abs(prod2.squeeze()[i]) < 1e-10

    @pytest.mark.parametrize(
            "w_ib_b, euler_1s, idx", [
                ( [20, 0, 0], [19.9, 0, 0], 0 ),
                ( [0, -20, 0], [0, -19.9, 0], 1 ),
                ( [0, 0, 5], [0, 0, 4.9], 2 ),
    ])
    def test_dqdt_signals(self, w_ib_b, euler_1s, idx):
        import scipy.integrate  as Int

        #  I: inertial frame
        #  b: body frame
        qI2b = kArrayNav( [0,0,0] ).euler2Q()

        # angular rotation between I and b:
        w_ib_b = kArrayNav( w_ib_b, hvector=False ).to_rad()

        def eqdiff(q,t,w_ib_b):
            """
            we will use scipy to call this function, and therefore q shall be a list().
            """
            qi2b = kArrayNav(q)
            dqdt = qi2b.dqdt( w_ib_b )
            return dqdt.to_list()

        T = np.linspace(0,1,20)
        y = Int.odeint(eqdiff, qI2b.to_list(), T, (w_ib_b,))[-1]

        # from the quaternions to euler, to compare:
        euler = kArrayNav(y, hvector=False).Q2euler().to_deg()

        if euler_1s[idx] > 0:
            assert euler_1s[idx] < euler[0][idx] < (euler_1s[idx] + 0.2)
        else:
            assert euler_1s[idx] > euler[0][idx] > (euler_1s[idx] - 0.2)

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
            q_inv = q.q_inv()

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

    def test_Mplus(self):
        for _ in range(50):
            euler_rad = kArrayNav( self.get_random_euler_rad() )
            qa2b = kArrayNav( euler_rad).euler2Q()

            euler_rad = kArrayNav( self.get_random_euler_rad() )
            qb2c = kArrayNav( euler_rad ).euler2Q()

            # ground truth:
            qa2c = qa2b.q_x_q(qb2c)

            # qa2c = qa2b.M+ x qb2c
            qtst = qa2b.to_Mplus() * qb2c

            # test:
            for i,j in zip(qa2c, qtst):
                assert abs(i-j) < 1e-10

    def test_Mminus(self):
        for _ in range(50):
            euler_rad = kArrayNav( self.get_random_euler_rad() )
            qa2b = kArrayNav( euler_rad).euler2Q()

            euler_rad = kArrayNav( self.get_random_euler_rad() )
            qb2c = kArrayNav( euler_rad ).euler2Q()

            # ground truth:
            qa2c = qa2b.q_x_q(qb2c)

            # qa2c = qb2c.M- x qa2b
            qtst = qb2c.to_Mminus() * qa2b

            # test:
            for i,j in zip(qa2c, qtst):
                assert abs(i-j) < 1e-10

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
