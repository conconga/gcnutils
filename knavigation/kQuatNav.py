#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kQuatNav.py
Beschreibung: Some basic operations with quaternions.
Autor: Luciano Auguto Kruk
Erstellt am: 06.09.2025
Version: 1.0.0
Lizenz:
GitHub:
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import numpy as np
from .kArray    import kArray
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kQuatNav:
    """
    The operations described here are based on Titterton & Weston, Strapdown Inertial Navigation Theory.
    The convention here is to have quaternions with the real part as the first element:

            q = R + i.U + j.V + k.W
    """

    def q_inv(self):
        """
        Calculates the inverse of a quaternion. This is similar to inverting a
        transformation matrix.
        """
        shape = self.shape
        tmp   = self.squeeze().tolist()
        ret   = [tmp[0]] + [-i for i in tmp[1:]]
        return self.__class__( np.asarray(ret).reshape(shape) / sum( [i**2 for i in ret] ))

    def q_conj(self):
        """
        Calculates the conjugated quaternion of the input.

        **Remark**: the first element is the real part of the quaternion here: (a + b.i + c.j + d.k)
        """
        #return self.__class__( np.hstack(([self[0]] + [-i for i in self[1:]])), hvector=False )
        return self.__class__( [-i if idx >= 1 else i for idx,i in enumerate(self) ], hvector=False )

    def q_norm(self):
        """
        Normalizes the given quaternion.
        """
        return self.__class__( self / np.sqrt(sum( [i**2 for i in self] )))

    def q_x_q(self, q2):
        """
        Navigation -- multiplies two quaternions
        (Titterton (3.55)-(3.56))

        **Remark** Here there are some hints on how to cascade quaternions for multiple
        transformations. See the respective test-units for details on the implementation.

            euler_a2b  -->  Ca2b | qa2b
            euler_b2c  -->  Cb2c | qb2c

            Ca2c = Cb2c * Ca2b
            qa2c = qa2b o qb2c  <== where 'o' is the product `q_x_q()`.

        output: np.array quaternion q3=q1.q2
        """

        q1 = self.squeeze()
        assert len(q1) == 4

        if hasattr(q2, "squeeze"):
            q2 = q2.squeeze()

        assert len(q2) == 4

        q3 = np.array([
            (q1[0]*q2[0]) - (q1[1]*q2[1]) - (q1[2]*q2[2]) - (q1[3]*q2[3]),
            (q1[1]*q2[0]) + (q1[0]*q2[1]) - (q1[3]*q2[2]) + (q1[2]*q2[3]),
            (q1[2]*q2[0]) + (q1[0]*q2[2]) + (q1[3]*q2[1]) - (q1[1]*q2[3]),
            (q1[3]*q2[0]) + (q1[0]*q2[3]) - (q1[2]*q2[1]) + (q1[1]*q2[2]),
        ])

        return self.__class__( q3, hvector=False )

    def _q_x_q(self, q2):
        """
        Same as q_x_q(), but calculate with another equation, more theoretical.

        q = [s v^T]^T

        """
        # assumption: q = [real, imag^T]^T
        q1 = self.squeeze()
        q2 = q2.squeeze()

        s1 = q1[0]
        v1 = q1[1:].view(np.ndarray)

        s2 = q2[0]
        v2 = q2[1:].view(np.ndarray)

        # q1 x q2:
        s3 = (s1*s2) - (v1.dot(v2))
        v3 = (s1*v2) + (s2*v1) + np.cross(v1,v2)

        return self.__class__( np.hstack( (s3,v3) ), hvector=False )

    def q_x_3d(self, vector):
        """
        Resolves a vector in another frame.
        This is similar to Ca2b x vector, but directly using quaternions.
        (Titterton (3.57))

        **Remark** Note that Titterton's equation matches his convention of
        creating quaternions using euler angles from frames a^ to b^, but
        performing transformations from b^ to a^.

        Here we provide an operation that uses a quaternion from a^ to b^ and
        a vector transformation also from a^ to b^. Hence:
        rb = q_a2b*  o  ra  o  q_a2b

        The vectors rb and ra are quaternion-form vectors (real-part = 0),
        but this is not important for inputs and output, which are 3D only.

        Inputs:
            Q4 (a + b.i + c.j + d.k)  :  quaternions
            vector                    :  3D vector

        Returns
            vector  :   3D vector resolved in another frame

        """

        # To be complete, we shall include the next line to invert the
        # quaternion and have the correct input for the transformation without
        # the C matrix.
        q = self.q_inv().squeeze() # <== first, use the inverse; see comments in the docstring
        v = vector.squeeze()

        # checks:
        assert len(q) == 4
        assert len(v) == 3

        # augmenting v to be an imaginary quaternion:
        v = np.hstack((0, v))

        # conjugate q:
        q_j = q.q_conj()

        # transformed vector (without the real part):
        ret = self.__class__(q).q_x_q(v).q_x_q(q_j)[1:] # <= to remove the real part

        return ret

    def dqdt(self, w):
        """
        The derivative of the quaternions is $\\dot{q} = 1/2 .B(w).q$
        This funtion returns $\\dot{q}$.

        Inputs:
            Q4 (a + b.i + c.j + d.k)  :  quaternions
            w_ib_b                    :  [rad/s] angular velocity
                
        Returns:
            d(q4)dt

        """

        assert len(self.squeeze()) == 4

        K      = 1e1
        cq     = kArray( [i for i in self], hvector=False )
        cq2    = kArray( [i*i for i in self], hvector=False )
        W      = kArray( w, hvector=True )
        epslon = 1.0 - sum(cq2.to_list())

        B = kArray( [
            [   0, -W[0][0], -W[0][1], -W[0][2]],
            [W[0][0],     0,  W[0][2], -W[0][1]],
            [W[0][1], -W[0][2],     0,  W[0][0]],
            [W[0][2],  W[0][1], -W[0][0],     0]
        ] )

        dq = (0.5 * B * cq) + (K*epslon*cq)
        return self.__class__( dq )

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
