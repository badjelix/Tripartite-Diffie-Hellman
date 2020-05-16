#!/bin/bash/python3

import math
from pynitefields import *

class Point:

    def __init__(self,x,y):
        self.x = x
        self.y = y
    def isInfinity(self):
        return False

class PointAtInfinity(Point):

    def __init__(self, prime):
        gf = GaloisField(prime)
        super().__init__(gf[0],gf[0])
    def isInfinity(self):
        return True

class EllipticCurve:

    def __init__(self,a,b):
        self.a = a
        self.b = b


#Elliptic curve operations

def negatePoint(p):
    return Point(p.x, p.y - 2 * p.y)

def addPoint(p, q, curve):
    if p.isInfinity():
        return q
    elif q.isInfinity():
        return p
    elif p.x == q.x and p.y == q.y - 2 * q.y:
        return PointAtInfinity(p.x.p)
    elif p.x == q.x and q.y == p.y:
        slope = (3 * pow(p.x, 2) + curve.a) / (2 * p.y)
    else:
        slope = (q.y - p.y) / (q.x - p.x)
    x = slope ** 2 - p.x - q.x
    y = slope * (p.x - x) - p.y
    return Point(x,y)

def negatePoint(p):
    return Point(p.x,p.y - 2 * p.y)

def getBinary(integer):
    return [int(n) for n in bin(integer)[2:]]

def doubleAndAdd(p, k, curve):
    binary = getBinary(k)
    q = PointAtInfinity()
    v = p
    i = len(binary) - 1
    while i >= 0:
        if binary[i] == 1:
            q = addPoint(q, v, curve)
        v = addPoint(v, v, curve)
        i -= 1
    return q


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
    return value.y - p.y + slope * (p.x - value.x)


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

gf = GaloisField(1009)
p = Point(gf[8],gf[703])
q = Point(gf[49],gf[20])
s = Point(gf[0],gf[0])

def WeilPairing(p, q, s, order, curve):
    a = Miller(p,order,addPoint(q,s,curve),curve)
    b = Miller(p,order,s,curve)
    c = Miller(q,order,addPoint(p, negatePoint(s), curve),curve)
    d = Miller(q,order,negatePoint(s),curve)
    return a * d / (b * c)

print(WeilPairing(p,q,s,7,EllipticCurve(gf[37],0))) 
print(WeilPairing(q,p,s,7,EllipticCurve(gf[37],0))) 

