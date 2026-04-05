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
        : input: kArray = [phi, theta, psi]
        : parameter : phi   [rad]
        : parameter : theta [rad]
        : parameter : psi   [rad]
        : return: kArray = Q4
        """
        array = self.squeeze()
        assert len(array) == 3
        phi   = array[0]
        theta = array[1]
        psi   = array[2]

        half_phi   = 0.5*phi
        half_theta = 0.5*theta
        half_psi   = 0.5*psi;

        return self.__class__( [
            (cos(half_phi)*cos(half_theta)*cos(half_psi)) + (sin(half_phi)*sin(half_theta)*sin(half_psi)),
            (sin(half_phi)*cos(half_theta)*cos(half_psi)) - (cos(half_phi)*sin(half_theta)*sin(half_psi)),
            (cos(half_phi)*sin(half_theta)*cos(half_psi)) + (sin(half_phi)*cos(half_theta)*sin(half_psi)),
            (cos(half_phi)*cos(half_theta)*sin(half_psi)) - (sin(half_phi)*sin(half_theta)*cos(half_psi))
        ], hvector=False );

    def Q2euler(self):
        """
        Navigation -- from Q to euler.
        : input: kArray = Q4
        : output   : phi   [rad]
        : output   : theta [rad]
        : output   : psi   [rad]
        : return: kArray( [phi, theta, psi] )
        """

        q     = self.squeeze()
        assert len(q) == 4

        phi   = atan2(2.0*((q[2]*q[3])+(q[0]*q[1])), (q[0]**2.0)-(q[1]**2.0)-(q[2]**2.0)+(q[3]**2.0));
        psi   = atan2(2.0*((q[1]*q[2])+(q[0]*q[3])), (q[0]**2.0)+(q[1]**2.0)-(q[2]**2.0)-(q[3]**2.0));
        try:
            theta = asin(2.0*((q[0]*q[2])-(q[1]*q[3])));
        except ValueError:
            print("ERR: norm(Q) = {:f}".format(np.sqrt(np.sum(q**2))))
            theta = 0;

        return self.__class__( (phi, theta, psi) )

    def Q2C(self):
        """
        Navigation -- from Q to C.

        If Q represents the transformation from 'a' to 'b', the matrix
        'C' represents 'Ca2b'.

        : input    : kArray = q
        : output   : C
        """

        q = self.squeeze()
        assert len(q) == 4

        C = kArray( np.empty((3,3)) )
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

    def C2euler(self):
        """
        Navigation -- from C to (phi,theta,psi)[rad]
        """

        C = self
        assert C.shape == (3,3)
        assert(C[2,2] != 0)
        assert(C[0,0] != 0)
        assert(C[0,2]>=-1 and C[0,2]<=1)

        phi   = np.arctan2(C[1,2], C[2,2])
        theta = np.arcsin(-C[0,2])
        psi   = np.arctan2(C[0,1], C[0,0])

        return self.__class__( (phi, theta, psi) )

    def euler2C(self):
        """
        : Convert euler [rad] to C matrix.
        : order = [phi, theta, psi] [rad]
        """

        a = self.squeeze()
        assert len(a) == 3

        if False:
            phi    = a[0]
            theta  = a[1]
            psi    = a[2]
            cphi   = cos(phi)
            ctheta = cos(theta)
            cpsi   = cos(psi)
            sphi   = sin(phi)
            stheta = sin(theta)
            spsi   = sin(psi)

            #C = self.__class__( np.empty((3,3)) )
            C = self.__class__.empty((3,3))
            C[0][0] = cpsi*ctheta
            C[0][1] = spsi*ctheta
            C[0][2] = -stheta
            C[1][0] = (-spsi*cphi) + (cpsi*stheta*sphi)
            C[1][1] = (cpsi*cphi) + (spsi*stheta*sphi)
            C[1][2] = ctheta * sphi
            C[2][0] = (spsi*sphi) + (cpsi*stheta*cphi)
            C[2][1] = (-cpsi*sphi) + (spsi*stheta*cphi)
            C[2][2] = ctheta * cphi

            #return C

        return self.__class__( a ).euler2Q().Q2C()

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

    def q1_x_q2(self, q2):
        """
        Navigation -- multiplies two quaternions

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
            (q1[0]*q2[0])-(q2[1]*q1[1])-(q2[2]*q1[2])-(q2[3]*q1[3]),
            (q2[0]*q1[1])+(q2[1]*q1[0])+(q2[2]*q1[3])-(q2[3]*q1[2]),
            (q2[0]*q1[2])+(q2[2]*q1[0])-(q2[1]*q1[3])+(q2[3]*q1[1]),
            (q2[0]*q1[3])+(q2[3]*q1[0])+(q2[1]*q1[2])-(q2[2]*q1[1])
        ])

        return self.__class__( q3 )

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
        See Titerton (3-52).
        : input  : [phi, theta, psi] [rad]
        : return : d[phi, theta, psi]/dt [rad/s]
        """

        a = self.squeeze()
        assert len(a) == 3
        phi    = a[0]
        theta  = a[1]
        psi    = a[2]

        if isinstance(w, list):
            b = w
        elif isinstance(w, self.__class__):
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
