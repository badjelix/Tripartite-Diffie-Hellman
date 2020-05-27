#!/usr/bin/python3

import math
import random
from util import *

class Participant:
    def __init__(self, name, curve, P, Q, order):
        """ Name from set {'A','B','C'} of the participant """
        self.name = name
        """ Dictionary to store the public values from different participants """
        self.publicKeys = {}
        """ Generate a private key: integer k s.t. 1 < k < order - 1 """
        self.privateKey = random.randint(1,order)
        """ Generate public values to broadcast to other 2 participants """
        # FIXME Points cant be points at infinity, right?
        self.publicKeys['P_' + name] = doubleAndAdd(P, self.privateKey, curve)
        self.publicKeys['Q_' + name] = doubleAndAdd(Q, self.privateKey, curve)


    """ Sends a tuple with the participant's public values """
    def sendPublicKeys(self):
        return (self.publicKeys['P_' + self.name], self.publicKeys['Q_' + self.name])

    """ Stores public values from another participant """
    def getPublicKeys(self, participant, keys):
        self.publicKeys['P_' + participant] = keys[0]
        self.publicKeys['Q_' + participant] = keys[1]
