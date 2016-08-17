import os

from npfmt import *
from Raise import *
from IsType import *


def Pause() : raw_input()

def Purge() : 
	try : os.system('purge')
	except : pass

