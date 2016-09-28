import numpy as np
import time



def RandomSeed() : 
	seed = int(('%.11f' % time.time()).split('.')[1][2:])
	np.random.seed(seed)
	return seed
