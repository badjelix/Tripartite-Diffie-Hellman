from participant import *
from galoisfield import *

if __name__ == "__main__":
    print("\nFIELD 31")
    gf31 = GaloisField(31)
    gf31[4].print()
    gf31[10].print()
    (gf31[4] * gf31[10]).print()
    print("AGORA OS A SERIO")
    el31_4 = GFieldElement(31,1,[4])
    el31_10 = GFieldElement(31,1,[10])
    el31_4.print()
    el31_10.print()
    (el31_4 * el31_10).print()
    print("\nFIELD 2^4")
    gf2 = GaloisField(2,4,[1, 1, 0, 0, 1])
    gf2[4].print()
    gf2[10].print()
    (gf2[4] * gf2[10]).print()
    print("AGORA OS A SERIO")
    el2_4 = GFieldElement(2,4,[1, 1, 0, 0],[1, 1, 0, 0, 1])
    el2_10 = GFieldElement(2,4,[1, 1, 1, 0],[1, 1, 0, 0, 1])
    el2_4.print()
    el2_10.print()
    (el2_4 * el2_10).print()

    """gf = GaloisField(17)
    ec = EllipticCurve(gf[0],gf[7],gf)
    generator = Point(gf[15],gf[13])

    A = Participant(ec, generator, 18)
    B = Participant(ec, generator, 18)
    C = Participant(ec, generator, 18)"""
