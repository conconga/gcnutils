#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kArrayNav.py
Beschreibung: Mathematical manipulation of vectors and matrices as needed by navigation modelling, and more.
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
from math     import sin, sqrt, cos, tan, atan, atan2, asin, pi
from mpmath   import sec
from .kArray  import kArray
from .kNavLib import kNavLib
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kNavTransformations(kNavLib):

    def to_deg(self):
        return self.__class__(self * 180. / pi)
        #return val * 180. / pi

    def to_rad(self):
        return self.__class__(self * pi / 180.)
        #return val * pi / 180.

    def euler2Q(self):
        """
        Navigation -- from euler to Q.
        (Titterton (3.65))

        **Remark** Farrell represents quaternions with its real element at the
            last position, after the imaginary components.

        Input:
            euler angles (kArray)  :  [rad]  [phi, theta, psi]

        Return
            (kArray)  :  Q4 (a + b.i + c.j + d.k)
        """
        array = self.squeeze()
        assert len(array) == 3

        half_phi, half_theta, half_psi = (0.5*i for i in array)

        return self.__class__( [
            (cos(half_phi)*cos(half_theta)*cos(half_psi)) + (sin(half_phi)*sin(half_theta)*sin(half_psi)),
            (sin(half_phi)*cos(half_theta)*cos(half_psi)) - (cos(half_phi)*sin(half_theta)*sin(half_psi)),
            (cos(half_phi)*sin(half_theta)*cos(half_psi)) + (sin(half_phi)*cos(half_theta)*sin(half_psi)),
            (cos(half_phi)*cos(half_theta)*sin(half_psi)) - (sin(half_phi)*sin(half_theta)*cos(half_psi))
        ], hvector=False );

    def Q2euler(self):
        """
        Navigation -- from Q to euler.
        (Titterton (3.66))

        Input:
            Q4 (a + b.i + c.j + d.k)  :  quaternions

        Return:
            (kArray) [rad] [phi, theta, psi]
        """

        q     = self.squeeze()
        assert len(q) == 4

        c11 = (q[0]**2.0) + (q[1]**2.0) - (q[2]**2.0) - (q[3]**2.0) # C[0,0]
        c21 = 2.0 * ((q[1]*q[2]) + (q[0]*q[3]))                     # C[0,1] <= transpose
        c31 = 2.0 * ((q[1]*q[3]) - (q[0]*q[2]))                     # C[0,2] <= transpose
        c32 = 2.0 * ((q[2]*q[3]) + (q[0]*q[1]))                     # C[1,2] <= transpose
        c33 = (q[0]**2.0) - (q[1]**2.0) - (q[2]**2.0) + (q[3]**2.0) # C[2,2]

        phi   = atan2(c32, c33)
        psi   = atan2(c21, c11)
        theta = asin(-c31)

        return self.__class__( (phi, theta, psi) )

    def Q2C(self):
        """
        Navigation -- from Q to C.
        (Titterton (3.59))

        **Remark**: See the remark in the docstring of `euler2C()`. The equation (3.63)
            shows that Titterton uses also here the transpose of the rotation
            from the static frame to the rotating frame.

        If Q represents the transformation from 'a' to 'b', the matrix
        'C' represents 'Ca2b'.

        Input:
            Q4 (a + b.i + c.j + d.k)  :  quaternions

        Return:
            C : transformation matrix

        """

        q = self.squeeze()
        assert len(q) == 4

        C = np.empty((3,3))
        C[0,0] = (q[0]**2.0) + (q[1]**2.0) - (q[2]**2.0) - (q[3]**2.0);
        C[0,1] = 2.0 * ((q[1]*q[2]) + (q[0]*q[3]));
        C[0,2] = 2.0 * ((q[1]*q[3]) - (q[0]*q[2]));

        C[1,0] = 2.0 * ((q[1]*q[2]) - (q[0]*q[3]));
        C[1,1] = (q[0]**2.0) - (q[1]**2.0) + (q[2]**2.0) - (q[3]**2.0);
        C[1,2] = 2.0 * ((q[2]*q[3]) + (q[0]*q[1]));

        C[2,0] = 2.0 * ((q[1]*q[3]) + (q[0]*q[2]));
        C[2,1] = 2.0 * ((q[2]*q[3]) - (q[0]*q[1]));
        C[2,2] = (q[0]**2.0) - (q[1]**2.0) - (q[2]**2.0) + (q[3]**2.0);

        return self.__class__( C )

    def _Q2C(self):
        """
        Same as Q2C(), but calculated using a theoretical equation:

        For a q = [s v^T], then

        C = (s^2 - v^T.v).I + 2.v.v^T -2.s.(v x)
        """
        q = self.squeeze()
        assert len(q) == 4

        s = q[0]
        v = self.__class__(q[1:], hvector=False)

        return (((s*s) - (v.T*v))*np.eye(3)) + (2 * v * v.T) - (2*s*v.to_skew())


    def C2euler(self):
        """
        Navigation -- from C to (phi,theta,psi)[rad]
        (Titterton (3.66))

        **Remark**: See the remark in the docstring of `euler2C()`. The approach here
            is described by Titterton (3.66), but with the transposed transformation
            matrix.

        Input:
            C : transformation matrix

        Return:
            (kArray) [rad] [phi, theta, psi]
        """

        c11 = self[0,0]
        c21 = self[0,1] # <= transpose
        c31 = self[0,2] # <= transpose
        c32 = self[1,2] # <= transpose
        c33 = self[2,2]

        phi   = atan2(c32, c33)
        psi   = atan2(c21, c11)
        theta = asin(-c31)

        return self.__class__( (phi, theta, psi) )

    def euler2C(self):
        """
        Convert euler [rad] to C matrix.
        (Titterton (3.63))

        **Remark**: the euler angles describe rotations of a frame (eg. body)
            over another (eg. navigation/tangent), and the transformation
            matrix is calculated by multiplying the individual rotations around
            Z-Y-X. Titterton does like this, but at (3.48) he transposes the
            transformation matrix and shows the transformation from the rotated
            frame to the static one. This method solves (3.47) instead, from the
            static to the rotated frame. See Farrell (2.31) for this result.

        Inputs:
            euler angles [phi, theta, psi] [rad]

        Return:
            C  :  Transformation Matrix
        """

        a = self.squeeze()
        assert len(a) == 3

        sphi, stheta, spsi = (sin(i) for i in a)
        cphi, ctheta, cpsi = (cos(i) for i in a)

        C = np.empty((3,3))
        C[0][0] = ctheta * cpsi
        C[0][1] = ctheta * spsi
        C[0][2] = -stheta

        C[1][0] = (sphi*stheta*cpsi) - (cphi*spsi)
        C[1][1] = (cphi*cpsi) + (sphi*stheta*spsi)
        C[1][2] = sphi*ctheta

        C[2][0] = (sphi*spsi) + (cphi*stheta*cpsi)
        C[2][1] = (cphi*stheta*spsi) - (sphi*cpsi)
        C[2][2] = cphi*ctheta

        return self.__class__( C )

    def C2Q(self):
        """
        Navigation -- from C to Q
        """
        return self.__class__( self.C2euler().euler2Q() )

    def Qinv(self):
        """
        Calculates the inverse of a quaternion. This is similar to inverting a transformation matrix.
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

    def q1_x_q2(self, q2):
        """
        Navigation -- multiplies two quaternions
        (Titterton (3.55)-(3.56))

        Let q1 represent C_a2b, and q2 represent C_b2c.
        The product C_a2c = C_b2c.C_a2b might be represented
        by q_a2c = q_b2c.q_a2b

        output: np.array quaternion q3=q1.q2
        """
        q1 = self.squeeze()
        assert len(q1) == 4

        if isinstance(q2, kArrayNav) or isinstance(kArray):
            q2 = q2.squeeze()
        assert len(q2) == 4

        q3 = np.array([
            (q1[0]*q2[0])-(q1[1]*q2[1])-(q1[2]*q2[2])-(q1[3]*q2[3]),
            (q1[1]*q2[0])+(q1[0]*q2[1])-(q1[3]*q2[2])+(q1[2]*q2[3]),
            (q1[2]*q2[0])+(q1[0]*q2[2])+(q1[3]*q2[1])-(q1[1]*q2[3]),
            (q1[3]*q2[0])+(q1[0]*q2[3])-(q1[2]*q2[1])+(q1[1]*q2[2])
        ])

        return self.__class__( q3, hvector=False )

    def q_x_3d(self, vector):
        """
        Resolves a vector in another frame.
        This is similar to Ca2b x vector, but directly using quaternions.
        (Titterton (3.57))

        The operation performed here is this:
        r2 = q . r1 . q*
        where r1 and r2 are vectors in quaternion-form (real part = 0), but
        with the inverse of `q`.

        With some algebra, the result above is the same as the operating the
        transformation backwards, or using the transpose of `Q2C()`. See the
        **remark** in the docstring of `Q2C()`.

        Inputs:
            Q4 (a + b.i + c.j + d.k)  :  quaternions
            vector                    :  3D vector

        Returns
            vector  :   3D vector resolved in another frame

        """

        q = self.q_inv().squeeze() # <== first, use the inverse
        v = vector.squeeze()

        # checks:
        assert len(q) == 4
        assert len(v) == 3

        # augmenting v to be an imaginary quaternion:
        v = np.hstack((0, v))

        # conjugate q:
        q_j = q.q_conj()

        # transformed vector (without the real part):
        ret = self.__class__(q).q1_x_q2(v).q1_x_q2(q_j)[1:] # <= to remove the real part

        return ret


    def ecef_llh2xyz(self):
        """
        : Convert from ECEF-geodetic to XYZe.
        : input   : llh  = (lat_rad, lon_rad, h_m)
        : output  : xzy_e [m]
        """

        geo = self.squeeze()
        assert len(geo) == 3
        lat = geo[0]
        lon = geo[1]
        h   = geo[2]

        s  = sin(lat);
        RN = self.earth_a / sqrt(1.0 - (self.earth_e2 * s * s));

        return self.__class__( [
            (RN + h) * cos(lat) * cos(lon),
            (RN + h) * cos(lat) * sin(lon),
            ((RN * (1.0 - self.earth_e2)) + h) * sin(lat)
        ], hvector=False )

    def ecef_xyz2llh(self):
        """
        : Convert from XYZe to ECEF-geodetic.
        : input  : xyz_e [m]
        : output : llh  = (lat_rad, lon_rad, h_m)
        """

        rect = self.squeeze()
        assert len(rect) == 3
        x = rect[0]
        y = rect[1]
        z = rect[2]

        p = sqrt((x * x) + (y * y));
        h     = 0;
        RN          = self.earth_a;
        for i in range(100): # timeout
            lasth   = h;

            # algorithm (Farrell/Barth p.28)
            s   = z / (((1.0 - self.earth_e2) * RN) + h);
            lat = atan((z + (self.earth_e2 * RN * s)) / p);
            RN  = self.earth_a / sqrt(1.0 - (self.earth_e2 * s * s));
            h   = (p / cos(lat)) - RN;

            # centimeter accuracy:
            if abs(lasth-h) < 0.01:
               break;

        lon = atan2(y, x);

        return self.__class__( [lat, lon, h] )

    def dqdt(self, w):
        """
        The derivative of the quaternions is $\\dot{q} = 1/2 .B(w).q$
        This funtion returns $\\dot{q}$.
        : input  : q4
        : output : d(q4)dt
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

    def dEulerDt(self, w):
        """
        Calculates the derivative vector of the euler angles.
        (Titterton (3-52))

        Inputs:
                euler angles : [rad]   [phi, thetha, psi]
            w : angular vel. : [rad/s] [wx, wy, wz]

        Return:
            d[phi, theta, psi]/dt [rad/s]
        """

        a = self.squeeze()
        assert len(a) == 3

        phi, theta, psi = a # [rad]

        if isinstance(w, list):
            b = w
        elif isinstance(w, (self.__class__, np.ndarray)):
            b = w.squeeze().tolist()
        else:
            raise(NameError("this type is still not supported"))

        aux = (b[1]*sin(phi)) + (b[2]*cos(phi))
        dphi = (aux * tan(theta)) + b[0]
        dtta = (b[1]*cos(phi)) - (b[2]*sin(phi))
        dpsi = aux * sec(theta)

        return self.__class__( [dphi, dtta, dpsi] )

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

class kArrayNav (kArray, kNavTransformations):
    def __init__(self, *args, **kargs):
        if len(args) > 0:
            super().__init__(*args, **kargs)

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
