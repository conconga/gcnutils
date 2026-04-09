#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kArrayLib.py
Beschreibung: Mathematical manipulation of vectors and matrices as needed by navigation modelling, and more.
Autor: Luciano Auguto Kruk
Erstellt am: 31.03.2026
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
from math import sqrt
from scipy.linalg import svd, qr
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kArrayLib:
    def to_skew(self):
        vtype = self._type(self)
        if vtype in [ self.TYPE_HORIZONTAL, self.TYPE_VERTICAL ]:
            val  = self.squeeze()
            assert len(val) == 3
            return self.__class__( [
                [0, -val[2], val[1]],
                [val[2], 0, -val[0]],
                [-val[1], val[0], 0] 
            ])
        else:
            raise(NameError("this is not a valid vector 3x1"))

    def X(self, y):
        """
        Calculates the cross-product between self and y.
        : input  : y = [3x1] vector
        : return : self X y
        """
        assert self.shape == (3,1)
        assert y.shape == (3,1)
        a = self.squeeze()
        b = y.squeeze()

        return self.__class__( [
            (a[1]*b[2]) - (a[2]*b[1]),
            (a[2]*b[0]) - (a[0]*b[2]),
            (a[0]*b[1]) - (a[1]*b[0])
        ], hvector=False )

    def apply(self, fn):
        """
        Applies a function to the array.
        : param : fn : function to call with the array as input
        : returns a new object with the result of 'fn'
        """

        return self.__class__( fn(self) )

    def norm(self):
        vtype = self._type(self)
        if vtype in [ self.TYPE_HORIZONTAL, self.TYPE_VERTICAL ]:
            return sqrt( sum( [i**2 for i in self.squeeze()] ))
        elif vtype == self.TYPE_SINGLEVALUE:
            return abs(float(self.squeeze()))
        else:
            raise(NameError(f"not prepared for type '{str(type(y))}' [self.__class__ = {self.__class__}]"))

    def inv(self):
        return np.linalg.inv(self).view(self.__class__)

    def to_orth(self, method="qr"):
        """
        Orthogonalizes the given matrix.

        The methods are:
           'qr'  :  decomposition QR
           'svd' :  nearest orthogonal matrix to 'self' in Frobenius norm (using SVD)
        """

        if method == "qr":
            Q, _ = qr(self, mode="economic")
            return self.__class__(Q)

        elif method == "svd":
            U, s, Vt = svd(self, full_matrices=False)
            return self.__class__(U @ Vt)

        else:
            raise NameError(f'method "{method}" is still not supported')

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
