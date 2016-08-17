
from Basic import *
import matplotlib as mpl
import matplotlib.pyplot as plt

from ShellCmd import *

name = ShellCmd('uname -a')[0][:11]
''' Use latex, may encounter error on node1 '''
if (name != 'Linux Node1') : mpl.rcParams['text.usetex'] = True
