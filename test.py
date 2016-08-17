#! /usr/bin/env python
from Basic import *
from Smooth import *
from Plot import *


a = np.random.random(10000)*2-1

b = 12.34
for i in xrange(10) : 
	t = np.linspace(0, np.random.random()*10*np.pi, a.size)
	b += np.sin(t)

n = (np.random.random(50)*a.size).astype(int)
c = np.random.random(n.size)*20

d = a + b
d[n] += c

e = Smooth(d, 0, 100, 100)

x = np.arange(a.size)
plt.plot(x, d)
plt.plot(x, b, 'c-', lw=3)
plt.plot(x, e, 'r-', lw=3)
plt.show()
