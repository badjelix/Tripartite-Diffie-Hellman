#!/usr/bin/python3

import math
import random
from util import *

class Participant:
    def __init__(self, name, curve, P, Q, order):
        """ Name from set {'A','B','C'} of the participant """
        self.name = name
        """ Dictionary to store the public values from different participants """
        self.publicValues = {}
        """ Generate a private key: integer k s.t. 1 < k < order - 1 """
        self.privateKey = random.randint(1,order)
        """ Generate public values to broadcast to other 2 participants """
        self.publicValues['P_' + name] = doubleAndAdd(P, self.privateKey, curve)
        self.publicValues['Q_' + name] = doubleAndAdd(Q, self.privateKey, curve)


        """s = doubleAndAdd(generator, random.randint(1,order), curve)

        while((s == self.P) or (s == self.Q) or (s == Point(0,0))):
            s.printPoint()
            s = doubleAndAdd(generator, random.randint(1,order), curve)

        weil = WeilPairing(self.P, self.Q, s, order, curve)"""
