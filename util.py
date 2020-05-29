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

    def testPoint(self,p):
        if p.y*p.y == p.x*p.x*p.x + self.a*p.x + self.b:
            return True
        else:
            return False


""" Elliptic curve operations """
def negatePoint(p):
    return Point(p.x, p.y - 2 * p.y)

def addPoint(p, q, curve):
    if p.isInfinity():
        return q
    elif q.isInfinity():
        return p
    elif p.x == q.x and p.y == q.y - 2 * q.y:
        return PointAtInfinity()
    elif p.x == q.x and q.y == p.y:
        slope = (3 * pow(p.x, 2) + curve.a) / (2 * p.y)
    else:
        slope = (q.y - p.y) / (q.x - p.x)
    x = slope ** 2 - p.x - q.x
    y = slope * (p.x - x) - p.y
    return Point(x,y)

def getBinary(integer):
    return [int(n) for n in bin(integer)[2:]]

def doubleAndAdd(p, n, curve):
    binary = getBinary(n)
    r = PointAtInfinity()
    q = p
    i = len(binary) - 1
    while i >= 0:
        if binary[i] == 1:
            r = addPoint(r, q, curve)
        q = addPoint(q, q, curve)
        i -= 1
    return r



# Miller Algorithm

def computeFunction(p, q, value, curve):
    if p.isInfinity() or q.isInfinity() or (p.x == q.x and p.y == q.y - 2 * q.y):
        if p.isInfinity():
            return value.x - q.x
        else:
            return value.x - p.x
    elif p.x == q.x and p.y == q.y:
        slope = (3 * pow(p.x, 2) + curve.a) / (2 * p.y)
    else:
        slope = (p.y - q.y) / (p.x - q.x)
    return (value.y - p.y + slope * (p.x - value.x)) / (value.x + p.x + q.x - pow(slope, 2))


def Miller(p, order, value, curve):
    res = 1
    v = p
    binary = getBinary(order)
    i = 1
    while i < len(binary):
        dv = addPoint(v, v, curve)
        res = res ** 2 * computeFunction(v, v, value, curve) #/ computeFunction(dv, negatePoint(dv), value, curve)
        v = dv
        if binary[i] == 1:
            vp = addPoint(v, p, curve)
            res = res * computeFunction(v, p, value, curve) #/ computeFunction(vp, negatePoint(vp), value, curve)
            v = vp
        i = i + 1
    return res

def WeilPairing(p, q, s, order, curve):
    a = Miller(p,order,addPoint(q,s,curve),curve)
    b = Miller(p,order,s,curve)
    c = Miller(q,order,addPoint(p, negatePoint(s), curve),curve)
    d = Miller(q,order,negatePoint(s),curve)

    return a * d / (b * c)

def TatePairing(p, q, order, curve, mod, n):
    a = Miller(p, order, q, curve)
    res = squareAndMultiply(a, ((mod ** n - 1) // order))
    return res

if __name__ == "__main__":

    prime = 47
    poly = [5,0,-4,0,1]
    ec = EllipticCurve(FieldElement([21,0,0,0],prime,4,poly),FieldElement([15,0,0,0],prime,4,poly))
    p = Point(FieldElement([45,0,0,0],prime,4,poly),FieldElement([23,0,0,0],prime,4,poly))
    q = Point(FieldElement([29,0,31,0],prime,4,poly),FieldElement([0,11,0,35],prime,4,poly))
    print(TatePairing(p,q,17,ec, prime, 4)) #expected: [39, 45, 43, 33]

    prime = 23
    poly = [1,0,1]
    ec = EllipticCurve(FieldElement([-1,0],prime,2,poly),FieldElement([0,0],prime,2,poly))
    p = Point(FieldElement([2,0],prime,2,poly),FieldElement([11,0],prime,2,poly))
    q = Point(FieldElement([21,0],prime,2,poly),FieldElement([0,12],prime,2,poly))
    s = Point(FieldElement([18,10],prime,2,poly),FieldElement([13,13],prime,2,poly))
    print(WeilPairing(p,q,s,3,ec)) #expected: [11, 15]

    p = Point(FieldElement(8,1009),FieldElement(703,1009))
    q = Point(FieldElement(49,1009),FieldElement(20,1009))
    s = Point(FieldElement(0,1009),FieldElement(0,1009))
    ec = EllipticCurve(FieldElement(37,1009),FieldElement(0,1009))
    print(WeilPairing(p,q,s,7,ec)) #expected: 105
