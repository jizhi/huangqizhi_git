#! /usr/bin/env python
import jizhipy as jp
from jizhipy.Plot import *
import numpy as np



sX = np.random.random()*90+10
a = np.random.randn(int(3e5)) * sX
dn = a.size/20


# jp.Smoothing() == Gaussian smoothing/blur/filter
per = np.arange(30)+3
times = np.arange(per.size)[::-1]+20


bstd1, bstd2 = [], []
for i in xrange(len(per)) : 
	print i+1, '/', len(per)
	yf, std = jp.GaussianFilter(per[i], times[i])

	xf = np.arange(yf.size)
	xf = xf - xf.mean()

	tf = (-200<=xf)*(xf<=200)
	xf = xf[tf]
	yf = yf[tf]

	xg = np.linspace(xf.min(), xf.max(), 3000)
	yg = jp.GaussianValue(xg, 0, std)

	plt.plot(xg, yg, 'b-', label='Gaussian(0, %.3f)' % std, lw=2)
	plt.plot(xf, yf, 'ro', label='per=%i, times=%i' % (per[i], times[i]), markersize=4)
	plt.legend(fontsize=13.5)
	plt.xlim(xf.min(), xf.max())
	plt.xlabel('x', size=16)
	plt.ylabel('Probability density', size=16)
	plt.title('Gaussian smoothing/blur/filter, (per=%i, times=%i), $\sigma=%.3f$' % (per[i], times[i], std))
	plt.savefig(('gaussian_smoothing_per-%i_times-%i_std-%.3f' % (per[i], times[i], std))+'.png')
	plt.close()


	b = jp.Smooth(a, 0, per[i], times[i])[dn:-dn]
	bstd1.append( b.std() )
	bstd2.append( a.std() / (std *2*np.pi**0.5)**0.5 )

bstd = np.array([bstd1, bstd2]).T
print bstd


