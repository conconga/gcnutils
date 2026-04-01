#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import sys
print( "**************************************" )
print(f"** __name__    = {__name__}")
print(f"** __package__ = {__package__}")
print(f"** sys.path[0] = {sys.path[0]}")
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

from knavigation import kArray
from unittest.mock import patch
import numpy as np
import pytest

if False:
    import pudb
    pudb.set_trace()


class TestClass_kArray:

    #----------------------------------------#
    @pytest.mark.singleshot
    def test_type_single_call(self):
        print("==== type() ====")

        # replacing the method "_type" by a mock:
        with patch.object(kArray, "_type") as mock_type:
            # configuring the return of the _type() method:
            mock_type.return_value = kArray.TYPE_ARRAY

            # construct the object with a matrix:
            a = kArray([[1,2],[3,4]])

            # asserts:
            assert mock_type.call_count == 1

    #----------------------------------------#
    @pytest.mark.parametrize(
            "content, vtype", [
                ( [1],           kArray.TYPE_SINGLEVALUE ),
                ( [1,2,3],       kArray.TYPE_HORIZONTAL ),
                ( [[1,2,3]],     kArray.TYPE_HORIZONTAL ),
                ( [[1]],         kArray.TYPE_SINGLEVALUE ),
                ( [[1,2]],       kArray.TYPE_HORIZONTAL ),
                ( [[1],[2]],     kArray.TYPE_VERTICAL ),
                ( [[1,2],[3,4]], kArray.TYPE_ARRAY ),
    ])
    def test_type_correct_classification(self, content, vtype):
        assert kArray._type( np.asarray( content ) ) == vtype

    #----------------------------------------#
    @pytest.mark.parametrize(
            "input_matrix, aux", [
                ( list(), 'just testing with a second argument' ),
                ( "string", "another test", ),
    ])
    def test_type_incorrect_classification(self, input_matrix, aux):
        print("input_matrix =", input_matrix)
        with pytest.raises(Exception): # <== ValueError?
            kArray._type( np.asarray( input_matrix ))

    #----------------------------------------#
    def test_vector_init(self):
        print("==== __init__() ====")

        print("kArray([1]) = ")
        print(kArray([1]))
        print("kArray([1,2,3]) = ")
        print(kArray([1,2,3]))
        print("kArray([1,2,3],hvector=True) = ")
        print(kArray([1,2,3],hvector=True))
        print("kArray([1,2,3],hvector=False) = ")
        print(kArray([1,2,3],hvector=False))
        print("kArray([[1,2,3]], hvector=False) = ")
        print(kArray([[1,2,3]], hvector=False))
        print("kArray([[1],[2]], hvector=True) = ")
        print(kArray([[1],[2]], hvector=True))

    def test_vector_transpose(self):
        print("==== transpose ====")
        a = kArray([1,2,3], hvector=True)
        print("a = {:f}".format(a))
        print("a.T = {:f}".format(a.T))
        print("a.T.T = {:f}".format(a.T.T))
        print("a.T.T.T = {:f}".format(a.T.T.T))
        print("a = {:f}".format(a))

    def test_vector_iter(self):
        print("==== iter ====")
        a = kArray( [1,2,3], hvector=False )

        # iter():
        for i,j in zip(a,[1,2,3]):
            assert i == j

    def test_vector_eq(self):
        print("==== eq ====")
        a = kArray([1,2,3])
        b = kArray([2,2,3])
        c = kArray([1,2])
        assert a == a
        assert a != b
        assert a != c
        a = kArray([8,9,10])
        print("a == a: {:s}".format((a==a).__str__()))
        print("a == b: {:s}".format((a==b).__str__()))

    def test_vector_sum(self):
        print("==== sum ====")
        a = kArray((1,2,3))
        b = kArray((4,5,6))
        assert a+b == kArray([5,7,9])
        assert a+(-1) == kArray([0,1,2])
        assert (-1)+a == kArray([0,1,2])
        print("a + b    = {:f}".format(a+b))
        print("a + (-1) = {:f}".format(a+(-1)))
        print("(-1) + a = {:f}".format((-1)+a))
        a += 7
        assert a == kArray([8,9,10])
        print("a += 7: {:f}".format(a))

    def test_vector_negative_signal(self):
        print("==== negative signal ====")
        a = kArray([8,9,10])
        print("-a = {:f}".format(-a))
        assert -a == kArray([-8,-9,-10])

    def test_vector_difference(self):
        print("==== difference ====")
        a = kArray((1,2,3))
        b = kArray((4,5,6))
        assert a-b == kArray([-3, -3, -3])
        assert 3.-a == kArray([2,1,0])
        print("a - b   = {:f}".format(a-b))
        print("3.0 - a = {:f}".format(3.0-a))
        a -= 3.0
        assert a == kArray([-2,-1,0])
        print("a -= 3.0: {:f}".format(a))

    def test_vector_multiplication(self):
        print("==== multiplication ====")
        a = kArray([1,2,3])
        b = kArray((4,5,6))
        try:
            print("a * b = {:f}".format(a*b))
        except:
            print("a * b =")
            print("    ^^^ERROR")

        assert a*b.T == 32
        assert 7.*a == kArray([7,14,21])
        print("a.T     = {:f}".format(a.T))
        print("7.0 * a = {:f}".format( 7.0 * a ))
        a *= 10
        assert a == kArray([10,20,30])
        print("a *= 10 : {:f}".format(a))

    #>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
    def test_matrix_init(self):
        print("==== __init__() ====")
        print("kArray( [[1,2,3],[4,5,6]] ) =")
        print(kArray( [[1,2,3],[4,5,6]] ))

    def test_matrix_transpose(self):
        print("==== transpose ====")
        a = kArray( [[1,2],[3,4]] )
        print("a = {:f}".format(a))
        print("a.T = {:f}".format(a.T))
        print("a.T.T = {:f}".format(a.T.T))
        print("a.T.T.T = {:f}".format(a.T.T.T))

        val = [1,3,2,4]
        for i in a.T:
            assert i == val.pop(0)

        assert a.T.T == a
        assert a.T.T.T == a.T

    def test_matrix_sum(self):
        print("==== sum ====")
        a = kArray( [[1,2,3], [4,5,6]] )
        b = kArray( [[0,1,2], [3,4,5]] )
        c = kArray( [[10,11],[12,13],[14,15]] )

        assert a+b == kArray( [[1,3,5], [7,9,11]] )
        assert c.T+a == kArray( [[11,14,17], [15,18,21]] )
        assert b+c.T == kArray( [[10,13,16], [14,17,20]] )

        assert a+(-1) == b
        assert (-1)+a == b

        a += 7
        assert a == kArray( [[8,9,10], [11,12,13]] )
        print("a += 7: {:f}".format(a))

    def test_matrix_negative_signal(self):
        print("==== negative signal ====")
        a = kArray( [[1,2,3], [4,5,6]] )
        assert -a == kArray( [[-1, -2, -3], [-4, -5, -6]] )
        print("-a = {:f}".format(-a))

    def test_matrix_difference(self):
        print("==== difference ====")
        a = kArray( [[1,2,3], [4,5,6]] )
        b = kArray( [[0,1,2], [3,4,5]] )
        c = kArray( [[10,11],[12,13],[14,15]] )

        assert a-b      == kArray( 1+np.zeros((2,3)) )
        assert a.T-c    == kArray( [[-9,-7],[-10,-8],[-11,-9]] )
        assert a-c.T    == kArray( [[-9,-10,-11],[-7,-8,-9]] )
        assert a.T-b.T  == kArray( np.ones((3,2)) )

        assert a-1 == kArray( [[0,1,2],[3,4,5]] )
        assert 1-a == kArray( [[0,-1,-2],[-3,-4,-5]] )

        a -= 3.0
        assert a == kArray( [[-2,-1,0],[1,2,3]] )
        print("a -= 3.0: {:f}".format(a))

    def test_matrix_multiplication(self):
        print("==== multiplication ====")
        a = kArray( [[1,2,3], [4,5,6]] ) # 2x3
        b = kArray( [[0,1,2], [3,4,5]] ) # 2x3
        c = kArray( [[10,11],[12,13],[14,15]] ) # 3x2
        d = kArray( [1,2,3], hvector=False ) # 3x1
        e = kArray( [2,3], hvector=True ) # 1x2

        assert a*b.T == kArray( [[8,26],[17,62]] )
        assert a.T*b == kArray( [[12,17,22],[15,22,29],[18,27,36]] )
        assert a*c == kArray( [[76, 82],[184,199]] )
        assert a*d == kArray( [[14],[32]] )
        assert e*a*d == 124.

        assert a*2.0 == kArray( [[2,4,6],[8,10,12]] )
        assert 2.0*b == kArray( [[0,2,4],[6,8,10]] )
        b *= -1.0
        assert b == kArray( [[0,-1,-2],[-3,-4,-5]] )

        try:
            ok = False
            print(a*a)
        except:
            ok = True

        if not ok:
            raise(NameError("it should not reach here"))

    def test_matrix_indexing(self):
        print("==== indexing ====")
        a = kArray( [[1,2,3], [4,5,6]] ) # 2x3
        assert a[1,1] == 5
        assert list(a[1]) == [4,5,6]

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
