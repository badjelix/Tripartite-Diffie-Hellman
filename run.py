from participant import *
from gfelement import *
from pythonSchoof.reduced_computation_schoof import *

if __name__ == "__main__":
    """print("\nFIELD 31")
    gf31 = GaloisField(31)
    gf31[4].print()
    gf31[10].print()
    (gf31[4] - gf31[10]).print()
    print("AGORA OS A SERIO")
    el31_4 = GFElement(31,1,[4])
    el31_10 = GFElement(31,1,[10])
    el31_4.print()
    el31_10.print()
    (el31_4 - el31_10).print()
    print("\nFIELD 2^4")
    gf2 = GaloisField(2,4,[1, 1, 0, 0, 1])
    gf2[4].print()
    gf2[10].print()
    (gf2[4] - gf2[10]).print()
    print("AGORA OS A SERIO")
    el2_4 = GFElement(2,4,[1, 1, 0, 0],[1, 1, 0, 0, 1])
    el2_10 = GFElement(2,4,[1, 1, 1, 0],[1, 1, 0, 0, 1])
    el2_4.print()
    el2_10.print()
    (el2_4 - el2_10).print()"""

    """ Generate the curve over a finite field"""
    field_prime = 29
    gf = GaloisField(field_prime)
    ec = EllipticCurve(gf[28],gf[1],gf)
    """ Get a generator of a large subgroup of the curve (h=2 or h=3) """
    curve_order = schoof.reduced_computation_schoof_algorithm(field_prime, 28, 1)
    print(curve_order)
    """ Generate 2 random independent points P and Q"""

    #A = Participant(ec, generator, 18)
    #B = Participant(ec, generator, 18)
    #C = Participant(ec, generator, 18)
