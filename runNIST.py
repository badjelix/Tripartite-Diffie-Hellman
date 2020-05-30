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
    p = 2**256 - 2**224 + 2**192 + 2**96 - 1
    k = 3
    q = p**k
    irreducible_poly = switchCoefs(gf_irreducible(k,p,ZZ))
    a = FieldElement([115792089210356248762697446949407573530086143415290314195533631308867097853948,0,0], p, k, irreducible_poly)
    b = FieldElement([41058363725152142129326129780047268409114441015993725554835256314039467401291,0,0], p, k, irreducible_poly)
    ec = EllipticCurve(a, b)

    """ Get 2 random independent points P and Q with the same order, by
        extending the field to k=3 so that p^k = 3mod4 and we can easily obtain
        square roots of elements of the field """
    # Just to make sure our assumption is correct
    if (q % 4) != 3:
        raise ValueError("q is not 3 mod 4")

    order = 115792089210356248762697446949407573529996955224135760342422259061068512044369

    Px = FieldElement([48439561293906451759052585252797914202762949526041747995844080717082404635286,0,0], p, k, irreducible_poly)
    Py = FieldElement([36134250956749795798585127919587881956611106672985015071877198253568414405109,0,0], p, k, irreducible_poly)
    P = Point(Px, Py)

    Qx = FieldElement([0,random.randint(1,order),0], p, k, irreducible_poly)
    Qy = sqrt3mod4(Qx*Qx*Qx+ec.a*Qx+ec.b,q)
    Q = Point(Qx, Qy)

    while(doubleAndAdd(Q, order, ec) != PointAtInfinity()):
        print("still no Q")
        Qx = FieldElement([0,random.randint(1,order),0], p, k, irreducible_poly)
        Qy = sqrt3mod4(Qx*Qx*Qx+ec.a*Qx+ec.b,q)
        Q = Point(Qx, Qy)

    print("JA TEMOS Q")

    S = Point(FieldElement([0,0,0], p, k, irreducible_poly), FieldElement([0,0,0], p, k, irreducible_poly))

    """ Create participants (each of them generates their private keys and public values) """
    print("A is generating his private and public keys")
    A = Participant('A', ec, P, Q, S, order)
    print("B is generating his private and public keys")
    B = Participant('B', ec, P, Q, S, order)
    print("C is generating his private and public keys")
    C = Participant('C', ec, P, Q, S, order)

    """ Broadcast phase: each participant sends their public values to the
        other two in a single message """
    # A sends public keys
    print("A is broadcasting his public values")
    B.getPublicKeys('A', A.sendPublicKeys())
    C.getPublicKeys('A', A.sendPublicKeys())
    # B sends public keys
    print("B is broadcasting his public values")
    A.getPublicKeys('B', B.sendPublicKeys())
    C.getPublicKeys('B', B.sendPublicKeys())
    # C sends public keys
    print("C is broadcasting his public values")
    A.getPublicKeys('C', C.sendPublicKeys())
    B.getPublicKeys('C', C.sendPublicKeys())

    """ Test generated key! """
    print("Participant A is going to broadcast a message to the other participants, type it here: (max 64 chars)")
    message = input()
    B.receiveMessage(A.sendMessage(message))
    C.receiveMessage(A.sendMessage(message))
