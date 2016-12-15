from scipy import interpolate
from npfmt import *
from Raise import *



def _Interp1d( xdata, ydata, xnew, kind='linear' ) : 
	'''
	1D interpolation. Note that all array must be 1D
	Outside the xdata, will use  Linear interpolation

	xdata:
		Must real

	ydata:
		Can be real and complex

	kind:
		'linear' or 'cubic'
	'''
	xdata, ydata, xnew = npfmt(xdata), npfmt(ydata), npfmt(xnew)
	if (len(xdata.shape)!=1 or len(ydata.shape)!=1 or len(xnew.shape)!=1) : Raise(Exception, 'xdata.shape='+str(xdata.shape)+', ydata.shape='+str(ydata.shape)+', xnew.shape='+str(xnew.shape)+' must all be 1D')
	#--------------------------------------------------
	xdata = xdata + 1j*ydata
	xdata = np.sort(xdata)
	ydata, xdata = xdata.imag, xdata.real  # xdata from min to max
#	ydata = np.arcsinh(ydata)  # smoother
	#--------------------------------------------------
	xnew = xnew + 1j*np.arange(xnew.size)
	xin = xnew[(xdata.min()<=xnew.real)*(xnew.real<=xdata.max())]
	xl  = xnew[(xdata.min()>xnew.real)]
	xr  = xnew[(xdata.max()<xnew.real)]
	#--------------------------------------------------
	funcin = interpolate.interp1d( xdata, ydata, kind=kind )
	yin = funcin(xin.real)
	#--------------------------------------------------
	for i in xrange(1, len(xdata)) : 
		if (xdata[i] != xdata[0]) : break
	yl = (ydata[i]-ydata[0])/(xdata[i]-xdata[0]) * (xl.real-xdata[0]) + ydata[0]
	#--------------------------------------------------
	for i in xrange(-2, -len(xdata), -1) : 
		if (xdata[i] != xdata[-1]) : break
	yr = (ydata[i]-ydata[-1])/(xdata[i]-xdata[-1]) * (xr.real-xdata[-1]) + ydata[-1]
	#--------------------------------------------------
	ynew = np.zeros(xnew.size)
	ynew[xin.imag.astype(int)] = yin
	ynew[xl.imag.astype(int)] = yl
	ynew[xr.imag.astype(int)] = yr
#	ynew = np.sinh(ynew)
	return ynew





def Interp1d( xdata, ydata, xnew, kind='linear' ) : 
	xdata, ydata, xnew = npfmt(xdata).flatten(), npfmt(ydata).flatten(), npfmt(xnew).flatten()
	if (ydata.dtype.name[:7] == 'complex') : 
		yr = _Interp1d(xdata, ydata.real, xnew, kind)
		yi = _Interp1d(xdata, ydata.imag, xnew, kind)
		ynew = yr + 1j*yi
	else : 
		ynew = _Interp1d(xdata, ydata, xnew, kind)
	return ynew
