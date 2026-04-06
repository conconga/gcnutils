#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kArrayNav.py
Beschreibung: Some standard equations based on WGS-84 model of the earth.
Autor: Luciano Auguto Kruk
Erstellt am: 08.10.2025
Version: 1.0.0
Lizenz: 
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
from math import sqrt, sin, cos, tan
import numpy as np
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kNavLib:
    # Earth Elliptic Model #
    earth_a  = 6378137.0; # [m]
    earth_b  = 6356752.3142; # [m]
    wie      = 1.0 * 7.2921151467e-5;
    earth_f  = (earth_a-earth_b)/earth_a;
    earth_e  = sqrt(earth_f*(2.0-earth_f));
    earth_e2 = (earth_e**2.0);

    def __init__(self, **kargs):
        super().__init__(**kargs)

    @classmethod
    def Rlambda(cls, lat_rad):
        """
        : parameter : lat_rad [rad] latitude
        : output    : R_lbd
        : (Farrell/Barth, eq 6-13)
        """
        return (cls.earth_a*(1.-cls.earth_e2)) / ((1.-(cls.earth_e2*(sin(lat_rad)**2)))**1.5);

    @classmethod
    def Rphi(cls, lat_rad):
        """
        : parameter : lat_rad [rad] latitude
        : output    : R_phi
        : (Farrell/Barth, eq 6-14)
        """
        return cls.earth_a / sqrt(1.-(cls.earth_e2*(sin(lat_rad)**2.0)));

    @classmethod
    def gravity(cls, lat_rad, h_m):
        """
        Using the geografic reference frame, this model allows:
            gravity_n = [0,0,output]
        (see Farrell, (6-141) and (6-142))

        : parameter : lat_rad [rad]  latitude
        : parameter : h_m     [m]    altitude above sea level
        : return    : gravity [m/s2] model of gravity for WGS-84
        """
        s2  = (sin(lat_rad))**2
        s22 = (sin(2.0*lat_rad))**2
        gamma_lat   = 9.780327 * ( 1. + (0.0053024*s2) - (0.0000058*s22) )
        gamma_lat_h = gamma_lat - ((3.0877e-6 - (0.0044e-6*s2) )*h_m) + (0.072e-12*(h_m**2))

        return gamma_lat_h

    @classmethod
    def dLat_dt(cls, vN, lat_rad, h_m):
        """
        Calculates the derivative of the geografic latitude.

        : parameter : vN      : [m/s] velocity-north
        : parameter : lat_rad : [rad] latitude
        : parameter : h_m     : [m]   altitude above sea level
        : return    : d(latitude)/dt : [rad/s]
        """

        return vN / (cls.Rlambda(lat_rad) + h_m)

    @classmethod
    def dLong_dt(cls, vE, lat_rad, h_m):
        """
        Calculates the derivative of the longitude.

        : parameter : vE      : [m/s] velocity-east
        : parameter : lat_rad : [rad] latitude
        : parameter : h_m     : [m]   altitude above sea level
        : return    : d(longitude)/dt : [rad/s]
        """

        return vE / (cos(lat_rad) * (cls.Rphi(lat_rad) + h_m))

    @classmethod
    def Re2n(cls, lat, lon):
        """
        Navigation -- calculates Re2n(lat,lon)
        The result (Re2n) does not change the content of the current object.

        : input    : lat   [rad]
        : input    : lon   [rad]
        : output   : Re2n
        """

        Re2n = np.empty((3,3));
        Re2n[0,0] = -sin(lat)*cos(lon);
        Re2n[0,1] = -sin(lat)*sin(lon);
        Re2n[0,2] = cos(lat);
        Re2n[1,0] = -sin(lon);
        Re2n[1,1] = cos(lon);
        Re2n[1,2] = 0;
        Re2n[2,0] = -cos(lat)*cos(lon);
        Re2n[2,1] = -cos(lat)*sin(lon);
        Re2n[2,2] = -sin(lat);

        return cls( Re2n )

    @classmethod
    def gravity_n(cls, lat_rad, h_m):
        """
        Calculates the local gravity vector in the geografic frame (n).

        : parameter : lat_rad [rad]  latitude
        : parameter : h_m     [m]    altitude above sea level
        : return    : vector with local gravity [3x1]
        """

        return cls( [0,0,cls.gravity(lat_rad, h_m)], hvector=False )

    @classmethod
    def dLLH_dt(cls, vN, vE, vD, lat_rad, h_m):
        """
        Calculates the vector with the derivatives of the
        latitude, longitude and altitude.

        : parameter : vN      : [m/s] velocity-north
        : parameter : vE      : [m/s] velocity-east
        : parameter : vD      : [m/s] velocity-down
        : parameter : lat_rad : [rad] latitude
        : parameter : h_m     : [m]   altitude above sea level
        """

        return cls( [
            cls.dLat_dt(vN, lat_rad, h_m),
            cls.dLong_dt(vE, lat_rad, h_m),
            -vD ], hvector=False )

    @classmethod
    def w_ie_n(cls, lat_rad):
        """
        Returns the angular velocity of the earth over the inertial frame, described at 'n'.

        : parameter : lat_rad : [rad] latitude
        """

        wie = cls.wie
        return cls( [
            wie * cos(lat_rad),
            0.0,
            -wie * sin(lat_rad)
        ], hvector=False )

    @classmethod
    def w_en_n(cls, dLat_dt, dLong_dt, lat_rad):
        """
        Calculates the angular velocity of the navigation frame over the earth frame, described at 'n'.

        : parameter : dLat_dt   : [rad/s] derivative of latitude
        : parameter : dLong_dt  : [rad/s] derivative of longitude
        : parameter : lat_rad   : [rad]   latitude
        """

        return cls( [
            dLong_dt * cos(lat_rad),
            - dLat_dt,
            - dLong_dt * sin(lat_rad)
        ], hvector=False )

    @classmethod
    def dWen_dt(cls, lat_rad, vN, vE, vD, h_m, vNp, vEp):
        """
        Calculates the first derivative of w_en_n, i.e d(w_en_n)/dt.
        (The expressions were obtained with sympy!)

        : parameter : lat_rad [rad]  latitude
        : parameter : vN      [m/s]  velocity north
        : parameter : vE      [m/s]  velocity east
        : parameter : vD      [m/s]  velocity down
        : parameter : h_m     [m]    height
        : parameter : vNp     [m/s2] d(vN)/dt
        : parameter : vEp     [m/s2] d(vE)/dt
        """

        dLat     = cls.dLat_dt(vN, lat_rad, h_m)
        slat     = sin(lat_rad)
        clat     = cos(lat_rad)
        e2s2l2   = cls.earth_e2 * slat**2.0

        wp_en_x =  ((cls.earth_a + sqrt(-e2s2l2 + 1.0)*h_m)*(-e2s2l2 + 1.0)*vEp - (1.0*cls.earth_a*cls.earth_e2*slat*clat*dLat + (-e2s2l2 + 1.0)**(3/2)*-vD)*vE)/((cls.earth_a + sqrt(-e2s2l2 + 1.0)*h_m)**2*sqrt(-e2s2l2 + 1.0))
        wp_en_y =  ((cls.earth_a*(cls.earth_e2 - 1.0) - (-e2s2l2 + 1.0)**1.5*h_m)*(-e2s2l2 + 1.0)**1.5*vNp - (-e2s2l2 + 1.0)**0.5*(1.5*cls.earth_a*cls.earth_e2*(cls.earth_e2 - 1.0)*sin(2*lat_rad)*dLat - (-e2s2l2 + 1.0)**2.5*-vD)*vN)/(cls.earth_a*(cls.earth_e2 - 1.0) - (-e2s2l2 + 1.0)**1.5*h_m)**2
        wp_en_z =  (-(cls.earth_a + sqrt(-e2s2l2 + 1.0)*h_m)*(-e2s2l2 + 1.0)*vE*slat**2*dLat - (cls.earth_a + sqrt(-e2s2l2 + 1.0)*h_m)*(-e2s2l2 + 1.0)*vE*clat**2*dLat - (cls.earth_a + sqrt(-e2s2l2 + 1.0)*h_m)*(-e2s2l2 + 1.0)*slat*clat*vEp + (1.0*cls.earth_a*cls.earth_e2*slat*clat*dLat + (-e2s2l2 + 1.0)**(3/2)*-vD)*vE*slat*clat)/((cls.earth_a + sqrt(-e2s2l2 + 1.0)*h_m)**2*sqrt(-e2s2l2 + 1.0)*clat**2)

        return cls( [ wp_en_x, wp_en_y, wp_en_z ], hvector=False )


    @classmethod
    def Jacobian_dwin_vNED(cls, lat_rad, h_m):
        """
        Calculates the Jacobian of w_in_n by vNED:
        (by sympy!)

            [ d(win_n_x)/dvN   d(win_n_x)/dvE   d(win_n_x)/dvD ]
        J = | d(win_n_y)/dvN   d(win_n_y)/dvE   d(win_n_y)/dvD |
            [ d(win_n_z)/dvN   d(win_n_z)/dvE   d(win_n_z)/dvD ]

        Input:
          lat_rad    :    [rad] latitude
          h          :    [m]   altitude

        """

        Rphi = cls.Rphi(lat_rad)
        Rlbd = cls.Rlambda(lat_rad)

        J = cls( [
            [     0,          1/(h_m + Rphi),          0   ],
            [-1/(h_m + Rlbd),         0,               0   ],
            [     0,    -tan(lat_rad)/(h_m + Rphi),    0   ],
        ])

        return J

    @classmethod
    def Jacobian_dwin_LH(cls, vN, vE, lat_rad, h_m):
        """
        Calculates the Jacobian of w_in_n by (lat/h):
        (by sympy!)

            [ d(win_n_x)/dlat    d(win_n_x)/dh ]
        J = | d(win_n_y)/dlat    d(win_n_y)/dh ]
            [ d(win_n_z)/dlat    d(win_n_z)/dh ]

        Input:
          vN      : [m/s] velocity-north
          vE      : [m/s] velocity-east
          lat_rad : [rad] latitude
          h_m     : [m]   altitude above sea level

        """
        slat  = sin(lat_rad)
        clat  = cos(lat_rad)
        s2lat = sin(2*lat_rad)
        c2lat = cos(2*lat_rad)
        slat2 = slat**2
        clat2 = clat**2
        lat   = lat_rad
        h     = h_m

        J_00 = -1.0*cls.earth_a*cls.earth_e2*vE*slat**1.0*clat/((cls.earth_a + h*sqrt(-cls.earth_e2*slat2 + 1.0))**2*sqrt(-cls.earth_e2*slat2 + 1.0)) - cls.wie*slat

        J_10 = -1.06066017177982*cls.earth_a*cls.earth_e2*vN*(cls.earth_e2 - 1.0)*(cls.earth_e2*c2lat - cls.earth_e2 + 2)**0.5*s2lat/(-cls.earth_a*cls.earth_e2 + cls.earth_a + h*(0.5*cls.earth_e2*c2lat - 0.5*cls.earth_e2 + 1)**1.5)**2

        J_20 = (1.0*cls.earth_a*cls.earth_e2*vE*slat2 - vE*(cls.earth_a + h*sqrt(-cls.earth_e2*slat2 + 1.0))*(-cls.earth_e2*slat2 + 1.0)/clat2 - cls.wie*(cls.earth_a + h*sqrt(-cls.earth_e2*slat2 + 1.0))**2*sqrt(-cls.earth_e2*slat2 + 1.0)*clat)/(
                (cls.earth_a + h*sqrt(-cls.earth_e2*slat2 + 1.0))**2*sqrt(-cls.earth_e2*slat2 + 1.0))

        J_01 = vE*(cls.earth_e2*slat2 - 1.0)/(cls.earth_a + h*sqrt(-cls.earth_e2*slat2 + 1.0))**2

        J_11 = vN*(-cls.earth_e2*slat2 + 1.0)**3.0/(cls.earth_a*(cls.earth_e2 - 1.0) - h*(-cls.earth_e2*slat2 + 1.0)**1.5)**2

        J_21 = vE*(-cls.earth_e2*slat2 + 1.0)*tan(lat)/(cls.earth_a + h*sqrt(-cls.earth_e2*slat2 + 1.0))**2

        J = cls( [
            [J_00, J_01],
            [J_10, J_11],
            [J_20, J_21],
        ])

        return J


#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
