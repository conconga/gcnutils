#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")

from knavigation import kArray
import numpy as np
from math import sqrt

class TestClass_kArrayLib:
    """
    The tests here are based on kArray, which is
    derived from kArrayLib.
    """

    def test_to_skew(self):
        print("==== to_skew() ====")
        a = kArray( [1,2,3], hvector=False )
        b = [0,-3,2,3,0,-1,-2,1,0]
        assert all([i==j for i,j in zip(a.to_skew(),b)])
        print(a.to_skew())

    def test_cross_product(self):
        print("==== cross product ====")
        a = kArray( [5,-2,4], hvector=True )
        b = kArray( [-1,2,3], hvector=True )

        try:
            ok = False
            c  = a.X(b)
        except:
            ok = True
        if not ok:
            raise(NameError("should not get here!"))

        c1 = a.T.X(b.T)
        c2 = np.cross(a.squeeze(), b.squeeze())
        c3 = a.to_skew() * b.T

        for i,j,k in zip(c1,c2,c3):
            assert abs(i-j) < 1e-10
            assert abs(i-k) < 1e-10

    def test_apply(self):
        print("==== apply() ====")
        a = kArray( [[2,0,0],[0,4,0],[0,0,10]] )
        b = a.apply(lambda x: np.linalg.inv(x))
        for i in range(3):
            assert b[i,i] == 1./a[i,i]


    def test_vector_norm(self):
        print("==== norm ====")
        a = kArray([1,1])
        assert abs(a.norm() - sqrt(2)) < 1e-10
        print("||a|| = {:f}".format(abs(a)))
        a = kArray(-32)
        assert abs(abs(a) - 32.0) < 1e-10

    def test_matrix_inv(self):
        print("==== inv() ====")
        M = [[1,2,3], [40,5,6], [7,8,9]]
        #M = np.random.randn(3,3)
        a = kArray(M)
        b = kArray(M)
        identity = a.inv() * b
        for i in range(3):
            for j in range(3):
                if i == j:
                    assert abs(identity[i,j] - 1.0) < 1e-10
                else:
                    assert abs(identity[i,j]) < 1e-10

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
