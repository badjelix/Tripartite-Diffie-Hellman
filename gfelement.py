#!/usr/bin/python3
import math

class GFElement():
    def __init__(self, p, n, coefs, irre_poly = []):
        """ Attributes of an element of a Galois Field: GF(p^n) """
        self.p = p
        self.n = n
        self.dim = int(pow(p, n))
        self.coefs = coefs
        self.irre_poly = irre_poly


    """ Addition """
    def __add__(self, el):
        # Make sure we're in the same field!
        if (self.p != el.p) or (self.n != el.n):
            print("Error, cannot add elements from different fields!")
            return None
        # Prime case
        if self.n == 1:
            return GFElement(self.p, self.n, [(self.coefs[0] + el.coefs[0]) % self.p])
        # Power of prime
        else:
            new_coefs = [(self.coefs[i] + el.coefs[i]) % self.p for i in range(0, self.n)]
            return GFElement(self.p, self.n, new_coefs, self.irre_poly)

    """ Add a field element to the left of this one """
    def __radd__(self, el):
        return self + el


    """ Subtraction """
    def __sub__(self, el):
        # Make sure we're in the same field!
        if (self.p != el.p) or (self.n != el.n):
            print("Error, cannot subtract elements from different fields!")
            return None
        # Prime case
        if self.n == 1:
            return GFElement(self.p, self.n, [(self.coefs[0] - el.coefs[0]) % self.p])
        # Power of prime
        else:
            new_coefs = [(self.coefs[i] - el.coefs[i]) % self.p for i in range(0, self.n)]
            return GFElement(self.p, self.n, new_coefs, self.irre_poly)


    """ Print out information about the element."""
    def print(self):
        print(self.coefs)













class GFElsement():
    def __init__(self, p, n, coefs, irre_poly = []):
        self.p = p
        self.n = n
        self.dim = int(pow(p, n))
        self.coefs = coefs
        self.irre_poly = irre_poly


    """ Multiplication """
    def __mul__(self, el):
        # Multiplication by a constant (must be on the right!)
        if isinstance(el, int):
            return GFElement(self.p, self.n, [(el * exp_coef) % self.p for exp_coef in self.coefs])
        # Multiplication by another FieldElement
        elif isinstance(el, GFElement):
            # Make sure we're in the same field!
            if (self.p != el.p) or (self.n != el.n):
                print("Error, cannot multiply elements from different fields!")
                return None
            # Prime case
            if self.n == 1:
                return GFElement(self.p, self.n, [(self.coefs[0] * el.coefs[0]) % self.p])
            # Power of prime case
            else:
                # Multiplying by 0, nothing to see here
                zeros = [0] * self.n
                if el.prim_power == zeros or self.prim_power == zeros:
                    return GFElement(self.p, self.n, zeros)
                else:
                    new_exp = self.prim_power + el.prim_power # New exponent
                    if new_exp > self.dim - 1:
                        new_exp = ((new_exp - 1) % (self.dim - 1)) + 1
                    new_coefs = [int(x) for x in self.field_list[new_exp].split(',')]
                    return GFElement(self.p, self.n, new_coefs, self.field_list, self.is_sdb, self.sdb_field_list)
        else:
            raise TypeError("Unsupported operator")


    """ Multiplication from the left """
    def __rmul__(self, el): # Implementing rmul so we can multiply on the left by integers
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

            # Prime
            if self.n == 1:
                if self.prim_power == 0:
                    print("Cannot divide by 0.")
                    return
            # Power of prime
            else:
                if self.field_list.index(self.str_rep) == 0:
                    print("Cannot divide by 0.")
                    return
            # Actually do the division
            return self * el.inv()


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
            return FieldElement(self.p, self.n, [int(math.pow(self.prim_power, exponent)) % self.p])
        # Power of prime case
        else:
            new_coefs = []
            # 0, and any element to the 0 is 0 by convention
            if self.prim_power == 0 or exponent == 0:
                new_coefs = [int(x) for x in self.field_list[0].split(',')]
            else:
                new_exp = self.prim_power * exponent
                if new_exp > self.dim - 1:
                    new_exp = ((new_exp - 1) % (self.dim - 1)) + 1
                new_coefs = [int(x) for x in self.field_list[new_exp].split(',')]
            return FieldElement(self.p, self.n, new_coefs, self.field_list, self.is_sdb, self.sdb_field_list)


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
        if self.coefs != el.coefs:
            return False
        if self.field_list != el.field_list:
            return False
        return True


    def inv(self):
        """ Compute the multiplicative inverse of a field element.

            Returns:
                The FieldElement that is the inverse of this one. All
                elements have a multiplicative inverse except for 0;
                if 0 is passed, prints error message and returns None.

            Note: The trace of an element can be invoked in two ways. One can
            do el.inv() or inv(el).
        """
        if self.n == 1: # Prime case - brute force :(
            if self.prim_power == 0:
                print("Error, 0 has no multiplicative inverse.")
                return

            for i in range(0, self.p):
                if (self.prim_power * i) % self.p == 1:
                    return FieldElement(self.p, self.n, [i])
        else: # Power of prime case
            if self.prim_power == 0:
                print("Error, 0 has no multiplicative inverse.")
                return
            # Last element is always 1 which is it's own inverse
            elif self.prim_power == self.dim - 1:
                return self
            # All other elements, find exponent which sums to dim - 1
            else:
                new_coefs = [int(x) for x in self.field_list[self.dim - self.prim_power - 1].split(',')]
                return FieldElement(self.p, self.n, new_coefs, self.field_list, self.is_sdb, self.sdb_field_list)


    def print(self):
        """ Print out information about the element."""
        print(self.coefs)
