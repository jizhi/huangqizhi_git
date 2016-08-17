import sys
#sys.path.append('/Users/huangqizhi/huangqizhi_git/jizhipy/')

import os
import numpy as np

from Raise import *
from npfmt import *
from IsType import *


def Pause() : raw_input()

def Purge() : 
	try : os.system('purge')
	except : pass

