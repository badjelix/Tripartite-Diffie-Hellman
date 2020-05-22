#!/usr/bin/python3

import math
import random
from util import *

class Participant:
    def __init__(self):
        self.privateKey = random.getrandbits(128)
        self.sharedKey = ""
