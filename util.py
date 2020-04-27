#!/bin/bash/python3

import math

class Point:

    def __init__(self,x,y):
        self.x = x
        self.y = y
    def isInfinity(self):
        return False

class PointAtInfinity(Point):

    def __init__(self):
        super().__init__(0,0)
    def isInfinity(self):
        return True

class EllipticCurve:

    def __init__(self,a,b):
        self.a = a
        self.b = b

def modulo(number, mod):
    z = complexDivision(number, mod) * mod
    return number - z

def complexDivision(number, divisor):
    q = number / divisor
    res = math.floor(q.real)
    if(isinstance(number, complex)):
        res = int(q.real)
        res += int(q.imag)*1j
    return res

# EECA

def xgcd(a, b):
    """return (g, x, y) such that a*x + b*y = g = gcd(a, b)"""
    x0, x1, y0, y1 = 0, 1, 1, 0
    if a.imag == 0 and a.real < 0:
        a = -a
    while a != 0:
        q, a, b = complexDivision(b, a), modulo(b, a), a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return b, x0, y0

def modularInverse(number, mod):
    g, x, _ = xgcd(number, mod)
    if g != 1:
        raise ValueError(str(g) + "Isnt mod a prime?")
    return modulo(x, mod)

#Elliptic curve operations

def negatePoint(p):
    return Point(p.x, -p.y)

def addPoint(p, q, curve, mod):
    if p.isInfinity():
        return q
    elif q.isInfinity():
        return p
    elif modulo(p.x,mod) == modulo(q.x,mod) and modulo(p.y,mod) == modulo(q.y,mod):
        return doublePoint(p, curve, mod)
    elif modulo(p.x,mod) == modulo(q.x,mod) and modulo(q.y,mod) == modulo(-p.y,mod):
        return PointAtInfinity()

    slope = modulo((p.y - q.y) * modularInverse(p.x - q.x, mod), mod)
    x = modulo((slope ** 2 - p.x - q.x), mod)
    y = modulo((slope * (p.x - x) - p.y), mod)
    return Point(x,y)

def doublePoint(p, curve, mod):
    if modulo(p.x,mod) == 0:
        return PointAtInfinity()
    slope = modulo((3 * p.x ** 2 + curve.a) * modularInverse(2 * p.y, mod), mod)
    x = modulo((slope ** 2 - 2 * p.x), mod)
    y = modulo((slope * (p.x - x) - p.y), mod)
    return Point(x,y)

def negatePoint(p):
    return Point(p.x,-p.y)

# Miller Algorithm

def getBinary(integer):
    return [int(n) for n in bin(integer)[2:]]

def computeFunction(p, q, value, curve, mod):
    if modulo(p.x,mod) == modulo(q.x,mod) and modulo(p.y,mod) == modulo(q.y,mod) and modulo(p.y,mod) != 0:
        return computeTangent(p, value, curve, mod)
    elif modulo(p.x,mod) == modulo(q.x,mod) and modulo(p.y,mod) == modulo(-q.y,mod):
        return computeVertical(p, value, mod)
    else:
        slope = modulo((p.y - q.y) * modularInverse(p.x - q.x, mod), mod)
        return modulo((value.y - p.y - slope * (value.x - p.x)) * modularInverse(value.x + p.x + q.x - slope **2, mod), mod)

def computeTangent(p, value, curve, mod):
    slope = modulo(((3 * p.x ** 2 + curve.a) * modularInverse(2 * p.y,mod)), mod)
    return modulo((-slope * value.x + value.y - p.y + slope * p.x) * modularInverse(value.x + 2 * p.x - slope**2, mod), mod)

def computeVertical(p, value, mod):
    return modulo((value.x - p.x), mod)

def Miller(p, order, value, curve, mod):
    res = 1
    v = p
    binary = getBinary(order)
    i = len(binary) - 1
    while i >= 0:
        res = modulo(res ** 2 * computeFunction(v, v, value, curve, mod), mod)
        v = doublePoint(p, curve, mod)
        if binary[i] == 1:
            res = modulo(res * computeFunction(v, p, value, curve, mod), mod)
            v = addPoint(v, p, curve, mod)
        i = i - 1
    return res

p = Point(8, 703)
q = Point(49,20)
s = Point(0,0)

def WeilPairing(p, q, s, order, curve, mod):
    a = Miller(p,7,addPoint(q,s,curve,mod),curve,mod)
    b = Miller(p,7,s,curve,mod)
    c = Miller(q,7,addPoint(p,negatePoint(s),curve,mod),curve,mod)
    d = Miller(q,7,negatePoint(s),curve,mod)
    

    e = Miller(p,7,q,curve,mod)
    f = Miller(q,7,p,curve,mod)
    print(modulo(e * modularInverse(f,mod),mod))
    return modulo(a * d * modularInverse(b * c, mod), mod)

print(WeilPairing(p,q,s,7,EllipticCurve(37,0),1009)) 
