from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 
from scipy.interpolate import interp1d
from npfmt import *



##################################################
##################################################
##################################################



def plt_usetex( tf ) : 
	'''
	Some system, 
		set usetex=True , use tex
		set usetex=False, don't use tex

	But also some system,
		set usetex=True , raise error
		set usetex=False, use tex!

	Therefore, set True/False dependenting on your system
	'''
	rcParams['text.usetex'] = bool(tf)



##################################################
##################################################
##################################################



def plt_legend() : 
	print ' 2     9    0/1'
	print ' 6    10    5/7'
	print ' 3     8     4'



##################################################
##################################################
##################################################


def plt_axisformat( axis, fmt='sci' ) : 
	'''
	axis: 'x' or 'y'

	fmt: 
		(1) '%d' or '%.2f' or '%.1e' and so on: 
			use plt.gca().xaxis.set_major_formatter(FormatStrFormatter(fmt))

		(2) 'sci'/'scientific' or 'plain'(turns off scientific notation): 
			use plt.ticklabel_format(axis=axis, style=fmt)
	'''
	if ('%' in fmt.lower()) : 
		if (axis.lower() == 'x') : plt.gca().xaxis.set_major_formatter(FormatStrFormatter(fmt))
		else : plt.gca().yaxis.set_major_formatter(FormatStrFormatter(fmt))
	else : plt.ticklabel_format(axis=axis, style=fmt, scilimits=(0,0))
	

##################################################
##################################################
##################################################



def plt_period( xlist, period, dmajor, fmt='%i' ) : 
	'''
	xlist:
		plt.plot(xlist, y)

	period:
		For o'clock, period=24 hour
		For phase, period=360 degree

	dmajor:
		Value of each major
	'''
	n = 1.*xlist[0] / dmajor
	if (n == int(n)) : n = int(n)
	else : n = 1 + int(n)
	x0 = n * dmajor
	x0 = np.arange(x0, xlist[-1]+dmajor/10., dmajor)
	x1 = list(x0 % period)
	for i in xrange(len(x1)) : x1[i] = (fmt % x1[i])
	plt.xticks(x0, x1)
	return [x0, x1]



##################################################
##################################################
##################################################



def plt_xaxes( xy, which='both', times=[], fmt='%i', fontsize=12, ticksize=8, direction='in', label2On=True, color='k', showminor=False, pcolor=False ) : 
	axes = plt.gca()
	times = npfmt(times)
	fontsize = npfmt(fontsize)
	if (which == 'major') : 
		if (times.size != 1) : raise Exception('xy="'+xy+'", which="major", times.size must =1')
		if (type(fmt) != str) : raise Exception('xy="'+xy+'", which="major", fmt must "%n.mf" or "%ni"')
		tickloc = MultipleLocator(times.flatten()[0])
		fmt = FormatStrFormatter(fmt)
		axes.xaxis.set_major_locator(tickloc)
		axes.xaxis.set_major_formatter(fmt)
		plt.tick_params(which='major', direction=direction, colors=color, size=ticksize)
	elif (which == 'minor') : 
		if (times.size != 1) : raise Exception('xy="'+xy+'", which="minor", times.size must =1')
		if (type(fmt) != str) : raise Exception('xy="'+xy+'", which="minor", fmt must "%n.mf" or "%ni"')
		tickloc = MultipleLocator(times[0])
		axes.xaxis.set_minor_locator(tickloc)
		plt.tick_params(which='minor', direction=direction, colors=color, size=ticksize/2)
		if (showminor) : 
			fmt = FormatStrFormatter(fmt)
			axes.xaxis.set_minor_formatter(fmt)
	elif (which == 'both') : 
		if (xy == 'x') : ns = '2'
		else : ns = '4'
		if (times.size != 2) : raise Exception('xy="'+xy+'", which="both", times.size must ='+ns)
		if (type(fmt) == str) : fmt = [fmt, fmt]
		elif (len(fmt) != 2) : raise Exception('xy="'+xy+'", which="both", fmt must ="%n.mf" or ["%n1.m1f", "%n2.m2f"]')
		tickloc0 = MultipleLocator(times.flatten()[0])
		tickloc1 = MultipleLocator(times.flatten()[1])
		fmt0 = FormatStrFormatter(fmt[0])
		fmt1 = FormatStrFormatter(fmt[1])
		axes.xaxis.set_major_locator(tickloc0)
		axes.xaxis.set_major_formatter(fmt0)
		axes.xaxis.set_minor_locator(tickloc1)
		if (showminor) : 
			axes.xaxis.set_minor_formatter(fmt1)
		plt.tick_params(which='major', direction=direction, colors=color, size=ticksize)
		plt.tick_params(which='minor', direction=direction, colors=color, size=ticksize/2)
	axesa = axes.xaxis.get_major_ticks()
	axesi = axes.xaxis.get_minor_ticks()
	for tick in axesa:
		tick.label1.set_fontsize(fontsize)
		if (label2On) : 
			if (pcolor == False) : tick.label2On = True
			tick.label2.set_fontsize(fontsize)
	for tick in axesi:
		if (showminor) : 
			tick.label1.set_fontsize(int(fontsize*2/3.))
			tick.label2.set_fontsize(int(fontsize*2/3.))
		if (label2On) : 
			if (pcolor == False) : tick.label2On = True



def plt_yaxes( xy, which='both', times=[], fmt='%i', fontsize=12, ticksize=8, direction='in', label2On=True, color='k', showminor=False, pcolor=False) : 
	axes = plt.gca()
	times = npfmt(times)
	if (len(fmt) == 1) : fmt = fmt[0]
	if (which == 'major') : 
		if (times.size != 1) : raise Exception('xy="'+xy+'", which="major", times.size must =1')
		if (type(fmt) != str) : raise Exception('xy="'+xy+'", which="major", fmt must "%n.mf" or "%ni"')
		tickloc = MultipleLocator(times[0])
		fmt = FormatStrFormatter(fmt)
		axes.yaxis.set_major_locator(tickloc)
		axes.yaxis.set_major_formatter(fmt)
		plt.tick_params(which='major', direction=direction, colors=color, size=ticksize)
	elif (which == 'minor') : 
		if (times.size != 1) : raise Exception('xy="'+xy+'", which="minor", times.size must =1')
		if (type(fmt) != str) : raise Exception('xy="'+xy+'", which="minor", fmt must "%n.mf" or "%ni"')
		tickloc = MultipleLocator(times[0])
		axes.yaxis.set_minor_locator(tickloc)
		if (showminor) : 
			fmt = FormatStrFormatter(fmt)
			axes.yaxis.set_minor_formatter(fmt)
		plt.tick_params(which='minor', direction=direction, colors=color, size=ticksize/2)
	elif (which == 'both') : 
		if (xy == 'x') : ns = '2'
		else : ns = '4'
		if (times.size != 2) : raise Exception('xy="'+xy+'", which="both", times.size must ='+ns)
		if (type(fmt) == str) : fmt = [fmt, fmt]
		elif (len(fmt) != 2) : raise Exception('xy="'+xy+'", which="both", fmt must ="%n.mf" or ["%n1.m1f", "%n2.m2f"]')
		tickloc0 = MultipleLocator(times.flatten()[0])
		tickloc1 = MultipleLocator(times.flatten()[1])
		fmt0 = FormatStrFormatter(fmt[0])
		fmt1 = FormatStrFormatter(fmt[1])
		axes.yaxis.set_major_locator(tickloc0)
		axes.yaxis.set_major_formatter(fmt0)
		axes.yaxis.set_minor_locator(tickloc1)
		if (showminor) : 
			axes.yaxis.set_minor_formatter(fmt1)
		plt.tick_params(which='major', direction=direction, colors=color, size=ticksize)
		plt.tick_params(which='minor', direction=direction, colors=color, size=ticksize/2)
	axesa = axes.yaxis.get_major_ticks()
	axesi = axes.yaxis.get_minor_ticks()
	for tick in axesa:
		tick.label1.set_fontsize(fontsize)
		if (label2On) : 
			if (pcolor == False) : tick.label2On = True
			tick.label2.set_fontsize(fontsize)
	for tick in axesi:
		if (showminor) : 
			tick.label1.set_fontsize(int(fontsize*2/3.))
			tick.label2.set_fontsize(int(fontsize*2/3.))
		if (label2On) : 
			if (pcolor == False) : tick.label2On = True



def plt_axes( xy='both', which='both', times=[], fmt='%i', fontsize=12, ticksize=8, direction='in', label2On=False, color='k', showminor=False, pcolor=False) : 
	'''
	plt_axes() is used to set the axis ticks so that they look beautiful.

	xy:
		select axis, 'x' or 'y' or 'both'.

	which:
		'major' or 'minor' or 'both'.

	times:
		times of this value will show the tick. Can be one value, [,], np.array([,]) or the data itself.
		1. times = onve value, list, np.array.shape=(2,)
			(1) which='major' or 'minor', times is one value.
			(2) which='both', times is list or np.array [[xmajor, xminor], [ymajor,yminor]].
		2. If times = [np.array], means it is 2D. It means this np.array is the data to be ploted. For plt.plot(), it can be xarray or yarray. For plt.pcolormesh(), it is the map. It set times.shape=2D, then the real times will be calculated by the code automatically.

	fmt:
		format of the tick.

	fontsize:
		Control the size of the font.
		if (xy != 'both') : fontsize.size ==1
		elif (xy == 'both') : fontsize = [xsize, ysize] or one value

	ticksize:
		Control the size of the ticks, default=8.
		if (xy != 'both') : ticksize.size ==1
		elif (xy == 'both') : ticksize = [xticksize, yticksize] or one value

	direction:
		direction of the ticks, 'in' or 'out'

	label2On:
		=True or =False. Show the ticks on top and right of the figure.

	showminor:
		=False or =True. Show the minor ticks text or not.

	color:
		color of the ticks, default black.

	pcolor:
		use to plt.pcolormesh() or not. If True, there must be an error bar on the right, then we can't set yaxis.label2On=True
	'''
	if (fmt is None) : fmt = '%i'
	if (xy == 'x') : 
		plt_xaxes( xy=xy, which=which, times=times, fmt=fmt, fontsize=fontsize, ticksize=ticksize, direction=direction, label2On=label2On, color=color, showminor=showminor, pcolor=pcolor)
	elif (xy == 'y') : 
		plt_yaxes( xy=xy, which=which, times=times, fmt=fmt, fontsize=fontsize, ticksize=ticksize, direction=direction, label2On=label2On, color=color, showminor=showminor, pcolor=pcolor)
	elif (xy == 'both') : 
		times = npfmt(times)
		if (len(times) != 2) : raise Exception('xy=="both", len(times) must =2')
		if (type(fmt) == str) : fmt = [fmt, fmt]
		elif (len(fmt) != 2) : raise Exception('xy="both", which="both", fmt must ="%n.mf" or len(fmt) must =2')
		fontsize = npfmt(fontsize)
		if (fontsize.size == 1) : fontsize = [fontsize[0], fontsize[0]]
		ticksize = npfmt(ticksize)
		if (ticksize.size == 1) : ticksize = [ticksize[0], ticksize[0]]
		plt_xaxes( xy=xy, which=which, times=times[0], fmt=fmt[0], fontsize=fontsize[0], ticksize=ticksize[0], direction=direction, label2On=label2On, color=color, showminor=showminor, pcolor=pcolor)
		plt_yaxes( xy=xy, which=which, times=times[1], fmt=fmt[1], fontsize=fontsize[1], ticksize=ticksize[1], direction=direction, label2On=label2On, color=color, showminor=showminor, pcolor=pcolor)



##################################################
##################################################
##################################################



def plt_color( N, r2b=False, k=1 ) : 
	'''
	Use for plt.plot(x, y, color=)

	N:
		number of color from red to blue

	r2b:
		color from red to blue or blue to red?

	k:
		r, g, b = k*(r, g, b)
		Make the color darker
	'''
	if (N == 8) : 
		return ['k', (0.42,0,1), 'b','c','g', 'y', (1,0.5,0), 'r']
	elif (N == 7) : 
		return ['k', (0.42,0,1), 'b', 'c', 'g', 'y', 'r']
	elif (N == 6) : 
		return ['k', 'b', 'c', 'g', 'y', 'r']
	elif (N == 5) : 
		return ['k', 'b', 'c', 'g', 'r']
	elif (N == 4) : 
		return ['k', 'b', 'g', 'r']
	elif (N == 3) : 
		return ['b', 'g', 'r']
	elif (N == 2) : 
		return ['b', 'r']
	elif (N == 1) : 
		return ['b']
	#--------------------------------------------------
	_jet_data_my = npfmt([
		(255,   0,   0),
		(255, 128,   0),
		(200, 200,   0),
		(  0, 200,   0),
		(  0, 250, 250),
		(  0,   0, 255),
		(160,   0, 255)]) /255.
	r, g, b = _jet_data_my.T
	#--------------------------------------------------
	r, g, b = npfmt([r,g,b])*k
	x = np.linspace(0, 1, r.size)
	fr = interp1d(x, r)
	fg = interp1d(x, g)
	fb = interp1d(x, b)
	x = np.linspace(0, 1, N)
	r, g, b = fr(x), fg(x), fb(x)
	color = []
	for i in range(N) : 
		color = color + [(r[i], g[i], b[i])]
	if (r2b is False) : color = color[::-1]
	return color



##################################################
##################################################
##################################################



