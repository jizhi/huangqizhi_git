from npfmt import *
import numpy as np



def ArrayGroup( array, step=1 ) : 
	'''
	Devide array into several groups basing on step

	array:
		Will first flatten() and Sort() (small to large)

	step:
		array[i]-array[i-1]<=step, we will consider array[i] and array[i-1] are in the same group

	return:
		[arraygroup, arrayrange]

		arrayrange[:,0] (include), arrayrange[:,1] (exclude)
	'''
	array, group = np.sort(npfmt(array).flatten()), [0]
	for t in range(1, len(array)) : 
		if (abs(array[t]-array[t-1]) <= step) : 
			group.append(group[-1])
		else : group.append(group[-1]+1)
	group = np.array(group)
	arraygroup, arrayrange = [], []
	for t in range(group.max()+1) : 
		arraygroup.append(array[group==t])
		arrayrange.append([arraygroup[-1][0],arraygroup[-1][-1]+1])
	return [arraygroup, npfmt(arrayrange)]


