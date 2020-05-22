from participant import *

if __name__ == "__main__":
    gf = GaloisField(17)
    p  = Point(gf[9],gf[15])
    ec = EllipticCurve(gf[0],gf[7])

    A = Participant()
    B = Participant()
    C = Participant()

    doubleAndAdd(p,A.privateKey,ec).printPoint()
