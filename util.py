#!/usr/bin/python3

import math
from pynitefields import *

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
        self.discriminant = -16 * (4 * a*a*a + 27 * b*b)
        if not self.isSmooth():
            raise Exception("The curve %s is not smooth!" % self)

    def isSmooth(self):
        return self.discriminant != self.gf[0]

    def testPoint(self, x, y):
        return y*y == x*x*x + self.a * x + self.b


#Elliptic curve operations

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
    if (p.x == q.x and p.y == q.y - 2 * q.y) or p.isInfinity() or q.isInfinity():
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
    i = len(binary) - 2
    while i >= 0:
        dv = addPoint(v, v, curve)
        res = res ** 2 * computeFunction(v, v, value, curve) #/ computeFunction(dv, negatePoint(dv), value, curve)
        v = dv
        if binary[i] == 1:
            vp = addPoint(v, p, curve)
            res = res * computeFunction(v, p, value, curve) #/ computeFunction(vp, negatePoint(vp), value, curve)
            v = vp
        i = i - 1
    return res

def WeilPairing(p, q, s, order, curve):
    a = Miller(p,order,addPoint(q,s,curve),curve)
    b = Miller(p,order,s,curve)
    c = Miller(q,order,addPoint(p, negatePoint(s), curve),curve)
    d = Miller(q,order,negatePoint(s),curve)

    print(a)
    print(b)
    print(c)
    print(d)
    return a * d / (b * c)


## Main and other functions

if __name__ == "__main__":
    gf = GaloisField(2,8,[1, 0, 1, 1, 1, 0, 0, 0, 1])
    p = Point(gf[9],gf[15])
    #q = Point(gf[49],gf[20])
    #s = Point(gf[0],gf[0])
    ec = EllipticCurve(gf[0],gf[7])

    doubleAndAdd(p,15,ec).printPoint()
