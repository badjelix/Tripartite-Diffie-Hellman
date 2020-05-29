from participant import *
from fieldelement import *
import random
from sympy.polys.galoistools import gf_irreducible
from sympy.polys.domains import ZZ

if __name__ == "__main__":

    """ Here we are trying to use a well known NIST curve. We get a generator P
        of the curve (with cofactor = 1) which we know the order. We must still
        generate a random Q linearly independent of P with the same order """

    """ Generate the curve over a finite field"""
    p = 2**448 - 2 ** 224 - 1
    k = 1
    r = 181709681073901722637330951972001133588410340171829515070372549795146003961539585716195755291692375963310293709091662304773755859649779

    while(((p**k-1) % r) != 0):
        print(k)
        k += 1
    print("K IS" + str(k))
