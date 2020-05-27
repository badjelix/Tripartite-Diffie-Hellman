#!/usr/bin/python

import math
import numpy


def xgcd(a, b):
    """return (g, x, y) such that a*x + b*y = g = gcd(a, b)"""
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        #print("a: " + str(a) + " b: " + str(b))
        q, a, b = b // a, b % a, a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return b, x0, y0

def inverse(number, mod):
    g, x, _ = xgcd(number, mod)
    if g != 1:
        raise ValueError(str(g) + "Isnt mod a prime?")
    return x

def polyxgcd(a, b, irre, mod):
    """return (g, x, y) such that a*x + b*y = g = gcd(a, b)"""
    if degree(a) == 0:
        b, x0, y0 = xgcd(a[0], mod)
        b = [b % mod] + [0] * (degree(irre) - 1)
        x0 = [x0 % mod] + [0] * (degree(irre) - 1)
        y0 = [y0 % mod] + [0] * (degree(irre) - 1)
    elif degree(b) == 0:
        b, x0, y0 = xgcd(b[0], mod)
        b = [b % mod] + [0] * (degree(irre) - 1)
        x0 = [x0 % mod] + [0] * (degree(irre) - 1)
        y0 = [y0 % mod] + [0] * (degree(irre) - 1)
    else:
        order = max(len(a), len(b))
        x0, x1, y0, y1 = [0] * order, [1] + [0] * (order - 1), [1] + [0] * (order - 1), [0] * order
        while a != [0] * len(a):
            q, r = polydiv(b, a, mod)
            a, b = r, a
            y0, y1 = y1, polysub(y0, polymul(q, y1, irre, mod), mod)
            x0, x1 = x1, polysub(x0, polymul(q, x1, irre, mod), mod)
    rest = b[0]
    b = b[:degree(irre)]
    if rest != 1:
        for i in range(len(b)):
            b[i] = b[i] * inverse(rest, mod) % mod
            x0[i] = x0[i] * inverse(rest, mod) % mod
            y0[i] = y0[i] * inverse(rest, mod) % mod
    return b, x0, y0

def polyinverse(a, irre, mod):
    g, x, _ = polyxgcd(a, irre, irre, mod)
    if g != [1] + [0] * (len(g) - 1):
        raise ValueError(str(g) + " Isnt mod a prime? " + str(x))
    return x[:degree(irre)]

def polymul(a, b, irre, mod):
    new_exp_coefs = [0] * (len(a) + len(b) + 1)
    for e in range(len(a)):
        for i in range(len(b)):
            new_exp_coefs[e + i] = (new_exp_coefs[e + i] + a[e] * b[i]) % mod
    _ , r = polydiv(new_exp_coefs, irre, mod)
    return r[:degree(irre)]

def polysub(a, b, mod):
    order = max(len(a), len(b))
    res = []
    for i in range(order):
        if len(a) <= i:
            res += [b[i] % mod]
        elif len(b) <= i:
            res += [a[i] % mod]
        else:
            res += [(a[i] - b[i]) % mod]
    return res

def degree(poly):
    res = 0
    for e in range(len(poly)):
        if poly[e] != 0:
            res = e
    return res

def polydiv(poly, div, mod):
    quotient = [0] * len(poly)
    while degree(poly) >= degree(div) and poly != [0] * len(poly):
        temp = [0] * (degree(poly) - degree(div)) + div
        new_poly = []
        inv = inverse(temp[degree(poly)], mod)
        quotient[degree(poly) - degree(div)] = (quotient[degree(poly) - degree(div)] + inv * poly[degree(poly)]) % mod
        for e in range(len(temp)):
            temp[e] = temp[e] * inv * poly[degree(poly)] % mod
        poly = polysub(poly, temp, mod)
    return quotient, poly


class FieldElement():
    """ Class for an element in a finite field.

        Args:
            p (int): The prime order of the field this element is in.
            n (int): The degree of the field extension for the field this
                     element is in.
            exp_coefs (list): The set of expansion coefficients of this element
                              in terms of some basis.

        Attributes:
            p (int): The prime order of the field this element is in.
            n (int): The degree of the field extension for the field this
                     element is in.
            dim (int): The dimension of the field, :math:`p^n`.
            exp_coefs (list): The set of expansion coefficients of this element
                              in terms of the polynomial basis.
    """

    def __init__(self, exp_coefs, p, n = 1, irre_poly = []):
        self.p = p
        self.n = n
        self.dim = int(pow(p, n))

        # Set the expansion coefficients.
        # If we're in a prime field, the basis is 1, and
        # the coefficient is just the value
        self.exp_coefs = exp_coefs
        self.irre_poly = irre_poly


    def __add__(self, el):
        """ Addition.

            Args:
                el (FieldElement): A FieldElement to add to this one.

            Returns:
                A FieldElement which is this element + el. For prime fields
                this is simply addition modulo :math:`p`, for power-of-prime
                fields we must add using the exp_coefs.
        """
        # Make sure we're in the same field!
        if (self.p != el.p) or (self.n != el.n):
            print("Error, cannot add elements from different fields!")
            return None

        # Prime case
        if self.n == 1:
            return FieldElement((self.exp_coefs + el.exp_coefs) % self.p, self.p)
        else: # Power of prime case
            # Coefficients simply add modulo p
            new_coefs = [(self.exp_coefs[i] + el.exp_coefs[i]) % self.p for i in range(0, self.n)]
            return FieldElement(new_coefs, self.p, self.n, self.irre_poly)


    def __radd__(self, el):
        """ Add a field element to the left of this one.

            Addition in finite fields is commutative so this works just like
            the normal add. This is implemented so we can use 'sum'
            over lists of FieldElement.
        """
        return self + el


    def __sub__(self, el):
        """ Addition.

            Args:
                el (FieldElement): A FieldElement to subtract from this one.

            Returns:
                A FieldElement which is this element - el. For prime fields
                this is simply subtraction modulo :math:`p`, for power-of-prime
                fields we must subtract using the exp_coefs.
        """
        # Make sure we're in the same field!
        if (self.p != el.p) or (self.n != el.n):
            print("Error, cannot subtract elements from different fields!")
            return None

        # Prime case
        if self.n == 1:
            return FieldElement((self.exp_coefs - el.exp_coefs) % self.p, self.p)
        else:  # Power of prime case
            # Coefficients subtract modulo p
            new_coefs = [(self.exp_coefs[i] - el.exp_coefs[i]) % self.p for i in range(0, self.n)]
            return FieldElement(new_coefs, self.p, self.n, self.irre_poly)


    def __mul__(self, el):
        """ Multiplication.

            Args:
                el (int or FieldElement): An element to multiply with this one.
                      Can also pass an integer value.

            Returns:
                This element * el. For prime fields, this amounts to simple
                multiplication modulo :math:`p`.
        """
        # Multiplication by a constant (must be on the right!)
        if isinstance(el, int):
            if self.n == 1:
                return FieldElement((el * self.exp_coefs) % self.p, self.p)
            else:
                return FieldElement([(el * exp_coef) % self.p for exp_coef in self.exp_coefs], self.p, self.n, self.irre_poly)

        # Multiplication by another FieldElement
        elif isinstance(el, FieldElement):
            # Make sure we're in the same field!
            if (self.p != el.p) or (self.n != el.n):
                print("Error, cannot multiply elements from different fields!")
                return None

            # Prime case
            if self.n == 1:
                return FieldElement((self.exp_coefs * el.exp_coefs) % self.p, self.p)
            # Power of prime case
            else:
                return FieldElement(polymul(self.exp_coefs, el.exp_coefs, self.irre_poly, self.p), self.p, self.n, self.irre_poly)

        else:
            raise TypeError("Unsupported operator")


    def __rmul__(self, el): # Implementing rmul so we can multiply on the left by integers
        """ Multiplication from the left. """
        return self * el


    def __truediv__(self, el):
        """ Division.

            In a Galois Field division is just multiplying by the inverse. By
            definition of a finite field, every element has a multiplicative
            inverse, except for 0.

            Args:
                An element to divide this one by.

            Returns:
                This element / el. Returns None if el = 0.
        """
        if isinstance(el, FieldElement):
            if (self.p != el.p) or (self.n != el.n):
                print("Error, cannot divide elements from different fields.")
            return self * el.inv()


    # Operations with assignment
    def __iadd__(self, el):
        """ Addition with assignment. """
        return self + el


    def __isub__(self, el):
        """ Subtraction with assignment. """
        return self - el


    def __imul__(self, el):
        """ Multiplication with assignment. """
        return self * el


    def __itruediv__(self, el):
        """ Division with assignment. """
        return self / el


    def __pow__(self, exponent):
        """ Exponentiation.

            Args:
                exponent (int): Something to exponentiate this element by.

            Returns:
                This element to the power of exponent. Just the normal power
                modulo p for primes. For power-of-primes, we define that the
                power of any element to 0 is the 0 element, and *not* 1.
        """
        # Prime case
        if self.n == 1:
            return FieldElement(int(self.exp_coefs ** exponent) % self.p, self.p)
        # Power of prime case
        else:
            poly = numpy.polynomial.polynomial.polypow(numpy.array(self.exp_coefs), exponent)
            poly = [i for i in poly]
            if len(poly) < len(self.irre_poly) - 1:
                poly += [0] * (len(self.irre_poly) - len(poly) - 1)
            _ , r = polydiv(poly, self.irre_poly, self.p)
            return FieldElement(r[:degree(self.irre_poly)], self.p, self.n, self.irre_poly)

    def __eq__(self, el):
        """ Test equality of two field elements.

            Args:
                el (FieldElement): An element to compare with.

            Returns:
                True if the field dimensions (:math:`p`, :math:`n`) are the
                same, the basis expansions are the same, and the list of
                field elements is the same. False otherwise.
        """
        if (self.p != el.p) or (self.n != el.n):
            return False
        if self.exp_coefs != el.exp_coefs:
            return False
        return True


    def __lt__(self, el):
        """ Implement a 'natural' ordering for field elements.

            For prime fields, this is simply the ordering of natural numbers.
            For power of primes, turn the coefficient lists into binary
            strings, and order them this way. Doing this to allow for
            Wigner functions to be plotted 'in order' in Balthasar.

            Args:
                el (FieldElement): An element to compare with.

            Returns:
                True if this element is 'less' by the conditions defined above.
                False otherwise.
        """
        if self.n == 1:
            if self.exp_coefs < el.exp_coefs:
                return True
            else:
                return False
        else:
            # If there is a sdb defined, use that, otherwise use exp_coefs
            this_exp_str = [str(x) for x in self.exp_coefs]
            that_exp_str = [str(x) for x in el.exp_coefs]
            if "".join(this_exp_str) < "".join(that_exp_str):
                return True
            else:
                return False


    def __repr__(self):
        """ Make the field element get printed in the command line."""
        return str(self.exp_coefs)


    def __hash__(self):
        """ Make hashable so we can use these guys as dictionary keys."""
        return hash(repr(self))


    def inv(self):
        """ Compute the multiplicative inverse of a field element.

            Returns:
                The FieldElement that is the inverse of this one. All
                elements have a multiplicative inverse except for 0;
                if 0 is passed, prints error message and returns None.

            Note: The trace of an element can be invoked in two ways. One can
            do el.inv() or inv(el).
        """

        if self.n == 1: #prime case
            if self.exp_coefs == 0:
                print("Error, 0 has no multiplicative inverse.")
                return

            return FieldElement(inverse(self.exp_coefs, self.p), self.p)

        else: # Power of prime case
            return FieldElement(polyinverse(self.exp_coefs, self.irre_poly, self.p), self.p, self.n, self.irre_poly)


    def print(self):
        """ Print out information about the element."""
        print(self.exp_coefs)

    def toString(self):
        """ Return string with information about the element."""
        return str(self.exp_coefs)
