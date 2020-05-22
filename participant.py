#!/usr/bin/python3

import math
import random
from util import *

class Participant:
    def __init__(self, curve, generator, curveorder):
        self.privateKey = random.randint(1,curveorder)
        self.P = doubleAndAdd(generator, random.randint(1,curveorder), curve)
        self.Q = doubleAndAdd(generator, random.randint(1,curveorder), curve)
        self.sharedKey = ""

        self.P.printPoint()
        self.Q.printPoint()

        WeilPairing(self.P, self.Q, s, order, curve)
