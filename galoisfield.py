#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# galoisfield.py: A class for elements in a finite field inspired in pynitefields
#
import math

class GFieldElement():
    """ Class for an element in a finite field.

        Args:
            p (int): The prime order of the field this element is in.
            n (int): The degree of the field extension for the field this
                     element is in.
            exp_coefs (list): The set of expansion coefficients of this element
                              in terms of some basis.
            irre_poly (list): Irreducible polynomial from the field this element belongs to.

        Attributes:
            p (int): The prime order of the field this element is in.
            n (int): The degree of the field extension for the field this
                     element is in.
            dim (int): The dimension of the field, :math:`p^n`.
            exp_coefs (list): The set of expansion coefficients of this element
                              in terms of the polynomial basis.
            irre_poly (list): Irreducible polynomial from the field this element belongs to.
    """

    def __init__(self, p, n, exp_coefs, irre_poly = []):
        self.p = p
        self.n = n
        self.dim = int(pow(p, n))
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
            return GFieldElement(self.p, self.n, [(self.exp_coefs[0] + el.exp_coefs[0]) % self.p])
        else: # Power of prime case
            # Coefficients simply add modulo p
            new_coefs = [(self.exp_coefs[i] + el.exp_coefs[i]) % self.p for i in range(0, self.n)]
            return GFieldElement(self.p, self.n, new_coefs, self.irre_poly)


    def __radd__(self, el):
        """ Add a field element to the left of this one.

            Addition in finite fields is commutative so this works just like
            the normal add. This is implemented so we can use 'sum'
            over lists of FieldElement.
        """
        return self + el


    def __sub__(self, el):
        """ Subtraction.

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
            return GFieldElement(self.p, self.n, [(self.exp_coefs[0] - el.exp_coefs[0]) % self.p])
        else:  # Power of prime case
            # Coefficients subtract modulo p
            new_coefs = [(self.exp_coefs[i] - el.exp_coefs[i]) % self.p for i in range(0, self.n)]
            return GFieldElement(self.p, self.n, new_coefs, self.irre_poly)


    def __mul__(self, el):
        """ Multiplication.

            Args:
                el (int or FieldElement): An element to multiply with this one.
                      Can also pass an integer value.

            Returns:
                This element * el. For prime fields, this amounts to simple
                multiplication modulo :math:`p`. For power of primes, this is
                where the ugly field_list comes in handy. We can compute the
                new power of the primitive element by adding together this one
                and the one from el; we then use field_list to find the
                corresponding FieldElement and return it.
        """
        # Multiplication by a constant (must be on the right!)
        if isinstance(el, int):
            return GFieldElement(self.p, self.n, [(el * exp_coef) % self.p for exp_coef in self.exp_coefs])

        # Multiplication by another FieldElement
        elif isinstance(el, GFieldElement):
            # Make sure we're in the same field!
            if (self.p != el.p) or (self.n != el.n):
                print("Error, cannot multiply elements from different fields!")
                return None

            # Prime case
            if self.n == 1:
                return GFieldElement(self.p, self.n, [(self.exp_coefs[0] * el.exp_coefs[0]) % self.p])
            # Power of prime case
            else:
                # I stored the whole list of field elements in each element for a reason...
                # Now we can multiply really easily

                # Multiplying by 0, nothing to see here
                if el.prim_power == 0 or self.prim_power == 0:
                    zeros = [0] * self.n
                    return GFieldElement(self.p, self.n, zeros)
                else:
                    new_exp = self.prim_power + el.prim_power # New exponent
                    # If the exponent calculated is outside the range of primitive element
                    # powers of the field, we need to wrap it around using the fact that
                    # the last field element is 1.
                    if new_exp > self.dim - 1:
                        new_exp = ((new_exp - 1) % (self.dim - 1)) + 1
                    new_exp_coefs = [int(x) for x in self.field_list[new_exp].split(',')]
                    return GFieldElement(self.p, self.n, new_exp_coefs, self.field_list, self.is_sdb, self.sdb_field_list)
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
        if self.exp_coefs != el.exp_coefs:
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
        print(self.exp_coefs)
