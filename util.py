#!/usr/bin/python3

import math
from pynitefields import *

""" Definition of Points and of the Elliptic Curve itself """
class Point:
    def __init__(self,x=0,y=0):
        self.x = x
        self.y = y

    def isInfinity(self):
        return False

    def printPoint(self):
        if(self.isInfinity()):
            print("Point at Infinity")
        else:
            print("X: ", end ="")
            self.x.print()
            print("Y: ", end ="")
            self.y.print()

class PointAtInfinity(Point):
    def __init__(self):
        super().__init__()

    def isInfinity(self):
        return True

class EllipticCurve:
    def __init__(self,a,b,gf):
        self.a = a
        self.b = b
        self.gf = gf
        #self.discriminant = -16 * (4 * a*a*a + 27 * b*b)
        #if not self.isSmooth():
        #    raise Exception("The curve %s is not smooth!" % self)

    #def isSmooth(self):
    #    return self.discriminant != self.gf[0]

    #def testPoint(self, x, y):
    #    return y*y == x*x*x + self.a * x + self.b


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


def squareAndMultiply(x, n):
    binary = getBinary(n)
    result = x / x
    i = len(binary) - 1
    while i >= 0:
        if binary[i] == 1:
            result *= x
        x *= x
        i -= 1
    return result


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
    ec = EllipticCurve(FieldElement(prime, 4, [21,0,0,0], irre_poly=poly) ,FieldElement(prime, 4, [15,0,0,0], irre_poly=poly))
    p = Point(FieldElement(prime,4,[45,0,0,0], irre_poly=poly),FieldElement(prime,4,[23,0,0,0],irre_poly=poly))
    q = Point(FieldElement(prime,4,[29,0,31,0], irre_poly=poly),FieldElement(prime,4,[0,11,0,35],irre_poly=poly))
    print(TatePairing(p,q,17,ec, prime, 4))

    prime = 23
    poly = [1,0,1]
    ec = EllipticCurve(FieldElement(prime, 2, [-1,0], irre_poly=poly) ,FieldElement(prime, 2, [0,0], irre_poly=poly))
    p = Point(FieldElement(prime,2,[2,0], irre_poly=poly),FieldElement(prime,2,[11,0],irre_poly=poly))
    q = Point(FieldElement(prime,2,[21,0], irre_poly=poly),FieldElement(prime,2,[0,12],irre_poly=poly))
    s = Point(FieldElement(prime,2,[18,10], irre_poly=poly),FieldElement(prime,2,[13,13],irre_poly=poly))
    print(WeilPairing(p,q,s,3,ec))

    gf = GaloisField(1009)
    p = Point(gf[8],gf[703])
    q = Point(gf[49],gf[20])
    s = Point(gf[0],gf[0])
    ec = EllipticCurve(gf[37],gf[0])
    print(WeilPairing(p,q,s,7,ec))

"""
## Main and other functions

def getBinary(integer):
    return [int(n) for n in bin(integer)[2:]]

if __name__ == "__main__":
    gf = GaloisField(2,8,[1, 0, 1, 1, 1, 0, 0, 0, 1])
    p = Point(gf[9],gf[15])
    #q = Point(gf[49],gf[20])
    #s = Point(gf[0],gf[0])
    ec = EllipticCurve(gf[0],gf[7])

    doubleAndAdd(p,15,ec).printPoint()
    """
