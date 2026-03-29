#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kArray.py
Beschreibung: Mathematical manipulation of vectors and matrices.
Autor: Luciano Auguto Kruk
Erstellt am: 06.09.2025
Version: 1.0.0
Lizenz: 
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import numpy as np
#import copy
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kArrayCommon:

    TYPE_ARRAY       = 0
    TYPE_VERTICAL    = 1
    TYPE_HORIZONTAL  = 2
    TYPE_SINGLEVALUE = 3

    @classmethod
    def _type(cls, val):
        """
            TYPE_ARRAY       = 2D array
            TYPE_VERTICAL    = vertical vector
            TYPE_HORIZONTAL  = horizontal vector
            TYPE_SINGLEVALUE = single value array
        """

        assert isinstance(val, np.ndarray)
        assert all( [i>0 for i in val.shape] )

        size = val.shape
        if len(size) == 1:
            if size[0] == 1:
                return cls.TYPE_SINGLEVALUE
            else:
                return cls.TYPE_HORIZONTAL
        elif len(size) == 2:
            if size == (1,1):
                return cls.TYPE_SINGLEVALUE
            if size[0] == 1:
                return cls.TYPE_HORIZONTAL
            elif size[1] == 1:
                return cls.TYPE_VERTICAL
            else:
                return cls.TYPE_ARRAY
        else:
            print("::error::")
            print(val)
            raise(NameError("I am not prepared for this array"))

    def _do_format_2D(self, fmt, R, C):
        txt = "[ ["
        for r in range(R):
            for c in range(C):
                txt += "{{:{:s}}}".format(fmt).format(self[r,c])
                if c < (C-1):
                    txt += ", "
            if r < (R-1):
                txt += "], ["
            else:
                txt += "]"
        txt += " ]"
        return txt

    def _do_format(self, fmt):
        size  = self.shape
        txt = self._do_format_2D(fmt, *size)
        return txt

    #( --- miscelaneous --- )#
    def __eq__(self, y):
        tol = 1e-10
        ret = True
        y   = self.__class__(y)

        if self.shape != y.shape:
            ret = False

        if not all( [ abs(i-j) <= max( [abs(i), abs(j)] )*tol for i,j in zip(self, y) ] ):
            ret = False

        return ret

    def __ne__(self, y):
        return not self.__eq__(y)

    def to_list(self):
        return self.reshape(-1).tolist()

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kArray (kArrayCommon, np.ndarray):

    def __new__(cls, *args, hvector=None):
        """
        When 'val' is given as a list, 'hvector' is used to indicate whether the
        vector to be created is horizontal (row, True) or vertical (column, False).

        If no positional value is provided, the object will create a matrix [1x1] with a single zero.
        If 'hvector' is None, then the selection is automatic when possible.
        """

        assert len(args) in [0,1]
        if len(args) == 0:
            # empty matrix:
            val = [0,]
        else:
            val = args[0]

        if isinstance(val, (list, tuple)):
            val = np.asarray( val )
        elif isinstance(val, (int, float)):
            val = np.asarray( [val] )
        else:
            # i guess the input is already an array
            pass

        # input type, vertical or horizontal or single?
        vtype = cls._type(val)
        assert 1 <= len(val.shape) <= 2

        if vtype == cls.TYPE_ARRAY:
            obj = np.asarray(args[0]).view(cls)

        elif vtype in [ cls.TYPE_HORIZONTAL, cls.TYPE_VERTICAL ]:
            if hvector == True:
                obj = val.squeeze().reshape(1,-1).view(cls)
            elif hvector == False:
                obj = val.squeeze().reshape(-1,1).view(cls)
            else: # hvector==None
                if vtype == cls.TYPE_HORIZONTAL:
                    obj = val.squeeze().reshape(1,-1).view(cls)
                else:
                    obj = val.squeeze().reshape(-1,1).view(cls)

        elif vtype == cls.TYPE_SINGLEVALUE:
            obj = val.squeeze().reshape(1,1).view(cls)

        else:
            print("::error::")
            print(val)
            print(vtype)
            raise(NameError("what is this?"))

        # casting to float to avoid issues like [1,2,3]+3.0 ==> CastingError
        obj = obj.astype(float)

        return obj

    def __repr__(self):
        return "{:s} |{:s}|".format(str(self.__class__.__name__), self._do_format("f"))

    def __init__(self, *args, **kargs):
        pass

    def __array_finalize__(self, obj) -> None:
        """
        Stattdessen wird array_finalize von NumPy selbst bei
        View-/Slicing-Operationen aufgerufen und sollte Attribute von obj
        kopieren oder initialisieren.

        Wird aufgerufen bei Slicing, uvm.
        """

        if obj is None: return

        # This attribute should be maintained!
        default_attributes = {"attr": 1}
        self.__dict__.update(default_attributes)

    def __format__(self, fmt):
        return self._do_format(fmt)

    # NumPy implements matrix multiplication other than we learn in the school.
    # Here I will bring it back.
    def __mul__(self, y):
        if isinstance(y, (int, float)):
            ret = super().__mul__(y)

        elif isinstance(y, (self.__class__)) or True: # <-- before failing, try this anyway (eg. kArrayNav)
            axb = self @ y
            ret = axb.view(self.__class__)

            # last check: if `ret` is a matrix [1x1], then return a float:
            if ret.shape in [(1), (1,1)]:
                ret = float(ret.squeeze())

        else:
            raise(NameError(f"not prepared for type '{str(type(y))}' [self.__class__ = {self.__class__}]"))

        return ret

    def __rmul__(self, y):
        # if y were a matrix, then __mul__() would be called by the matrix.
        return self * y
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

    #( --- iter --- )#
    def __iter__(self):
        temp = self.reshape(-1).tolist()
        for i in temp:
            yield i

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
