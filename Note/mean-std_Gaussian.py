#! /usr/bin/env python
import numpy as np
import jizhipy as jp

'''
Gaussian noises:
	a: a.mean()=0, a.std()=astd
	b: b.mean()=0, b.std()=bstd  (== or != astd)

Auto-correlation: 
	chi-squared distribution with k=1
	A = <a a*>
		A.mean() =           astd**2
		A.std()  = sqrt{2} * astd**2 = sqrt{2} * A.mean()

Cross-correlation:
	
	C = <a b*>
		C.mean() = 0
		C.std()  = astd * bstd = sqrt{Astd * Bstd / 2}


a.mean()=amean, a.std()=astd
A.mean() = amean**2 + astd**2
'''


# Two Gaussian Noises
a = np.random.randn(1000000) * 3.2
b = np.random.randn(a.size)  * 3.2


# Auto-correlation
dn = 1000
Atmp, A = a*a, []
for i in xrange(Atmp.size/dn) : 
	A.append(Atmp[i*dn:(i+1)*dn])
A = np.array(A)

Btmp, B = b*b, []
for i in xrange(Btmp.size/dn) : 
	B.append(Btmp[i*dn:(i+1)*dn])
B = np.array(B)


amean = ('%13.10f' % a.mean())
astd  = ('%13.10f' % a.std())
astd2 = ('%13.10f' % a.std()**2)
astd2sqrt2 = ('%13.10f' % (2**0.5*a.std()**2))
Amean = ('%13.10f' % A.mean())
Astd  = ('%13.10f' % A.std())
print 'a.mean()   =', amean, '   a.std()            =', astd
print 'a.std()**2 =', astd2, '   sqrt{2}*a.std()**2 =', astd2sqrt2
print 'A.mean()   =', Amean, '   A.std()            =', Astd
print
bmean = ('%13.10f' % b.mean())
bstd  = ('%13.10f' % b.std())
bstd2 = ('%13.10f' % b.std()**2)
bstd2sqrt2 = ('%13.10f' % (2**0.5*b.std()**2))
Bmean = ('%13.10f' % B.mean())
Bstd  = ('%13.10f' % B.std())
print 'b.mean()   =', bmean, '   b.std()            =', bstd
print 'b.std()**2 =', bstd2, '   sqrt{2}*b.std()**2 =', bstd2sqrt2
print 'B.mean()   =', Bmean, '   B.std()            =', Bstd
print


# Cross-correlation
Ctmp, C = a*b, []
for i in xrange(Ctmp.size/dn) : 
	C.append(Ctmp[i*dn:(i+1)*dn])
C = np.array(C)


Cmean = ('%13.10f' % C.mean())
Cstd  = ('%13.10f' % C.std())
abstd = ('%13.10f' % (a.std()*b.std()))
ABstd2sqrt = ('%13.10f' % (A.std()*B.std()/2)**0.5)
print '                              a.std()*b.std()    =', abstd
print '                         sqrt{A.std()*B.std()/2} =', ABstd2sqrt
print 'C.mean()   =', Cmean, '   C.std()            =', Cstd
print
print 'A=<a a*>     B=<b b*>     C=<a b*>'


