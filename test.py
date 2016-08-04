#! /usr/bin/env python
from Basic import *

from Sort import *

#a = (np.random.random((2,3,4))*100).astype(int)
#np.save('a', a)
a = np.load('a.npy')
print a
print
print '----------------'
print

b = Sort(a, '[1,:,2]')
print b
print
print '----------------'
print
print a[:,0,:]
print a[:,1,:]
print a[:,2,:]
print
print b[:,0,:]
print b[:,1,:]
print b[:,2,:]
