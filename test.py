#! /usr/bin/env python
import jizhipy as jp
import numpy as np
from jizhipy.Plot import *

color = plt_color(16)

a = np.arange(len(color))
b = a*0+1
for i in xrange(len(color)) : 
	plt.plot([a[i]], [b[i]], color=color[i], ls='', marker='o', markersize=10)
plt.savefig('color.png')




#def Print( array, precision=6, suppress=True ) : 
#	'''
#	Format the printing of np.array
#	suppress:
#		=False, print 1.23e+4
#		=True,  print 12340.
#	'''
#	array = np.array(array)
#	precision = int(round(precision))
#	if (not suppress) : 
#		fmt = '%i.%ie' % (8+precision, precision)
#		formatter = {'all': lambda x: format(x, fmt)}
#	#	formatter = {'all': lambda x: (('%'+fmt) % x)}
#	else : formatter = None
#	default = np.get_printoptions()
#	np.set_printoptions(precision=precision, suppress=suppress, formatter=formatter)
#	print array
#	np.set_printoptions(**default)
#
#
#
#array = np.array([-12.34, -0.123, -0.000034, 0, 0.0345, 0.00012, 0.00001, 0.000009, 56.78, -1234.56]).reshape(5, 2)
#
#
#Print(array, suppress=True)
#print
#
#valid, idx = jp.SciNot(array)
#print valid
#print idx
#
#print(array)
#print
#
#
##import numpy as np
##
##
##def Print( a, precision=6, suppress=True ) : 
##	'''
##	Format the printing of np.array
##	suppress:
##		=False, print 1.23e+4
##		=True,  print 12340.
##	'''
##	precision = int(round(precision))
##	if (not suppress) : formatter = {'all': lambda x: format(x, '.'+str(precision)+'e')}
##	else : formatter = None
##	default = np.get_printoptions()
##	np.set_printoptions(precision=precision, suppress=suppress, formatter=formatter)
##	print a
##	np.set_printoptions(**default)
##
##
##
##Print(array)
##
##
##
##shape = array.shape
##array = array.flatten()
##
### Get the sign
##sign = np.sign(array)  # sign(0)=0
##
### Convert to abs
##array = abs(array)
##
### because sign(0)=0, convert 0 to 1
##array[array==0] = 1
##
##nlarge, nsmall = (array>=1), (array<1)  # bool, not int
##
### Use log10 to get the power index
##idxlarge = np.log10(array[nlarge]).astype(int)
##
##scalesmall = int(round(np.log10(array[nsmall].min())))-2
##array[nsmall] /= 10.**scalesmall
##idxsmall = np.log10(array[nsmall]).astype(int) + scalesmall
##array[nsmall] *= 10.**scalesmall
##
##idx = np.zeros(array.size, int)
##idx[nlarge], idx[nsmall] = idxlarge, idxsmall
##
##valid = sign * (array / 10.**idx)
##
##valid, idx = valid.reshape(shape), idx.reshape(shape)
##
##Print(valid)
##Print(idx)
##
##
