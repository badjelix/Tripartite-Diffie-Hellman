from participant import *

if __name__ == "__main__":
    gf = GaloisField(2,4,[1, 1, 0, 0, 1])
    ele13 = FieldElement(2,4,[1, 0, 1, 1],irre_poly=[1, 1, 0, 0, 1])
    ele14 = FieldElement(2,4,[1, 0, 0, 1],irre_poly=[1, 1, 0, 0, 1])
    print(gf[14] / gf[13])
    print(ele14 / ele13)
    #A = Participant()
    #B = Participant()
    #C = Participant()

    #doubleAndAdd(p,A.privateKey,ec).printPoint()
