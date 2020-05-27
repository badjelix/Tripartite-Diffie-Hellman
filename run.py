from participant import *
from fieldelement import *

if __name__ == "__main__":

    """ Generate the curve over a finite field"""
    p = 482677778157700435350444108563600470389539607291135742953085077414483299007817968457323051999107203153032937333023591271636050696817523671646492380723773419011
    k = 2
    irreducible_poly = [1,0,1]
    ec = EllipticCurve(FieldElement([1,0], p, k, irreducible_poly), FieldElement([0,0], p, k, irreducible_poly))

    """ Get 2 random independent points P and Q with the same order (these
        points are the generators) """
    P = Point(FieldElement([44190300200219570605979955052143576952357255515115686851170191818316842095486907625480884395317616863401927551006066189692708095924815897927498508535823262371,0], p, k, irreducible_poly), FieldElement([26090947680860922395540330613428690525406329616428470738073031338841260885477380307130420220342204765301865163480203757570223664606235381540801075563801118751,0], p, k, irreducible_poly))
    Q = Point(FieldElement([417418390151798179157327683814659014460849518350508436411447781417311430237331232958577456865429161040089806217226455983348248260335272068783343983410685645620, 0], p, k, irreducible_poly), FieldElement([0,85984079438328066829535503806402848425113755688042614534609435398882015068450504353865472815063531531657210019063972911218641810155964304683033635085838106425], p, k, irreducible_poly))
    order = 593917583375891588584754753148372137203682206097
    S = Point(FieldElement([0,0], p, k, irreducible_poly), FieldElement([0,0], p, k, irreducible_poly))

    """ Create participants (each of them generates their private keys and public values) """
    A = Participant('A', ec, P, Q, order)
    B = Participant('B', ec, P, Q, order)
    C = Participant('C', ec, P, Q, order)

    """ Broadcast phase: each participant sends their public values to the
        other two in a single message """
    # A sends public keys
    B.getPublicKeys('A', A.sendPublicKeys())
    C.getPublicKeys('A', A.sendPublicKeys())
    # B sends public keys
    A.getPublicKeys('B', B.sendPublicKeys())
    C.getPublicKeys('B', B.sendPublicKeys())
    # C sends public keys
    A.getPublicKeys('C', C.sendPublicKeys())
    B.getPublicKeys('C', C.sendPublicKeys())
