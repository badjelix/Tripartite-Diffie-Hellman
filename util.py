#!/usr/bin/python3

import math
from fieldelement import *

""" Definition of Points and of the Elliptic Curve itself """
class Point:
    def __init__(self,x=0,y=0):
        self.x = x
        self.y = y

    def isInfinity(self):
        return False

    def toString(self):
        if(self.isInfinity()):
            return "Point at Infinity"
        else:
            return "X: " + self.x.toString() + "\nY: " + self.y.toString()

class PointAtInfinity(Point):
    def __init__(self):
        super().__init__()

    def isInfinity(self):
        return True

class EllipticCurve:
    def __init__(self,a,b):
        self.a = a
        self.b = b

    def testPoint(self,P):
        if P.y*P.y == P.x*P.x*P.x + self.a*P.x + self.b:
            return True
        else:
            return False


""" Elliptic curve operations """

# inverse of a point (x,y) is (x,-y)
def negatePoint(P):
    return Point(P.x, -P.y)

# compute P + Q = R
def addPoint(P, Q, curve):
    if P.isInfinity():
        return Q
    elif Q.isInfinity():
        return P
    elif P.x == Q.x and P.y == -Q.y:
        return PointAtInfinity()
    elif P.x == Q.x and P.y == Q.y:
        slope = (3 * pow(P.x, 2) + curve.a) / (2 * P.y)
    else:
        slope = (Q.y - P.y) / (Q.x - P.x)
    x = slope ** 2 - P.x - Q.x
    y = slope * (P.x - x) - P.y
    return Point(x,y)

# binary representation of an integer (ex.: 4 -> [1,0,0])
def getBinary(integer):
    return [int(n) for n in bin(integer)[2:]]

# compute [n]P
def doubleAndAdd(P, n, curve):
    binary = getBinary(n)
    R = PointAtInfinity()
    i = len(binary) - 1
    while i >= 0:
        if binary[i] == 1:
            R = addPoint(R, P, curve)
        P = addPoint(P, P, curve)
        i -= 1
    return R



""" Miller Algorithm """

# compute the function f(value) passing through P and Q
def computeFunction(P, Q, value, curve):
    if P.isInfinity() or Q.isInfinity() or (P.x == Q.x and P.y == -Q.y):
        if P.isInfinity():
            return value.x - Q.x
        else:
            return value.x - P.x
    elif P.x == Q.x and P.y == Q.y:
        slope = (3 * pow(P.x, 2) + curve.a) / (2 * P.y)
    else:
        slope = (P.y - Q.y) / (P.x - Q.x)
    return (value.y - P.y + slope * (P.x - value.x)) / (value.x + P.x + Q.x - pow(slope, 2))

# compute fP(value), with divisor order(P) - order(O)
def Miller(P, order, value, curve):
    res = 1
    V = P
    binary = getBinary(order)
    i = 1
    while i < len(binary):
        dV = addPoint(V, V, curve)
        res = res ** 2 * computeFunction(V, V, value, curve)
        V = dV
        if binary[i] == 1:
            VP = addPoint(V, P, curve)
            res = res * computeFunction(V, P, value, curve)
            V = VP
        i = i + 1
    return res

# Compute the Weil pairing between P and Q, with [order]P = [order]Q = O
# The pairing is independent from the choice of the point S, as long as S != {O,P,-Q,P-Q}
#
# the pairing is computed as fP(Q+S) * fQ(-S) / (fP(S) * fQ(P-S))
def WeilPairing(P, Q, S, order, curve):
    a = Miller(P,order,addPoint(Q,S,curve),curve)
    b = Miller(P,order,S,curve)
    c = Miller(Q,order,addPoint(P, negatePoint(S), curve),curve)
    d = Miller(Q,order,negatePoint(S),curve)

    return a * d / (b * c)

# Compute the Tate pairing between P and Q, with [order]P = O
#
# the pairing is computed as fP(Q)^((p^n - 1) / 2)
def TatePairing(P, Q, order, curve):
    a = Miller(P, order, Q, curve)
    res = squareAndMultiply(a, ((P.x.p ** P.x.n - 1) // order))
    return res

if __name__ == "__main__":

    prime = 47
    poly = [5,0,-4,0,1]
    ec = EllipticCurve(FieldElement([21,0,0,0],prime,4,poly),FieldElement([15,0,0,0],prime,4,poly))
    P = Point(FieldElement([45,0,0,0],prime,4,poly),FieldElement([23,0,0,0],prime,4,poly))
    Q = Point(FieldElement([29,0,31,0],prime,4,poly),FieldElement([0,11,0,35],prime,4,poly))
    print(TatePairing(P,Q,17,ec)) #expected: [39, 45, 43, 33]

    prime = 23
    poly = [1,0,1]
    ec = EllipticCurve(FieldElement([-1,0],prime,2,poly),FieldElement([0,0],prime,2,poly))
    P = Point(FieldElement([2,0],prime,2,poly),FieldElement([11,0],prime,2,poly))
    Q = Point(FieldElement([21,0],prime,2,poly),FieldElement([0,12],prime,2,poly))
    S = Point(FieldElement([18,10],prime,2,poly),FieldElement([13,13],prime,2,poly))
    print(WeilPairing(P,Q,S,3,ec)) #expected: [11, 15]

    P = Point(FieldElement(8,1009),FieldElement(703,1009))
    Q = Point(FieldElement(49,1009),FieldElement(20,1009))
    S = Point(FieldElement(0,1009),FieldElement(0,1009))
    ec = EllipticCurve(FieldElement(37,1009),FieldElement(0,1009))
    print(WeilPairing(P,Q,S,7,ec)) #expected: 105
