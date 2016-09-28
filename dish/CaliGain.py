import os
import numpy as np
import jizhipy as jp
from jizhipy.Plot import *
import scipy.signal as spsn
##################################################



class CaliGain( object ) : 


	def __init__( self, antarray=None, masking=None, Nprocess=None, verbose=True ) : 
		self.Nprocess = Nprocess
		self.verbose = verbose
		if (antarray is not None) : self.antarray = antarray
		if (masking  is not None) : self.masking  = masking
		self.outdir = jp.Outdir((None,'file'), (0,'file'))


	def Window( self, size ) : 
		''' How many time points to calculate one std? '''
		N = int(round(self.antarray.vis.shape[0] / size))
		self.nsplit = np.linspace(0, self.antarray.vis.shape[0], N+1).round().astype(int)



	def Gainnu( self, Nline=20, legendsize=10, legendloc=1 ) : 
		''' Calculate g(nu) and plot '''
		outdir = self.outdir + 'Gainnu/'
		if (not os.path.exists(outdir)) : os.makedirs(outdir)
		if (self.verbose) : 
			print 'CaliGain.Gainnu: start @', jp.Time(1)
			progressbar = jp.ProgressBar('    completed:', len(self.nsplit)-1+self.antarray.visorder.size)
		# Read data
		vistype = self.antarray.vistype[:-1]
		if (vistype != 'auto') : 
			print '    Error: CaliGaint.Gainnu() should use "auto-correlation", but now AntArray.vistype="'+self.antarray.vistype+'"'
			if (self.verbose) : print 'CaliGain.Gainnu:  end  @', jp.Time(1)+'\n'
			exit()
		#--------------------------------------------------
		vis = self.antarray.vis[:,:,self.antarray.visorder].real
		if (len(vis.shape) == 2) : vis = vis[:,:,None]  # 3D
		# Replace masks
		vis[self.masking.mask] = self.masking.maskvalue
		vis = np.ma.MaskedArray(vis, self.masking.mask)  #@#@
		# Calculate
		self.gainnu = np.zeros((self.nsplit.size-1,)+vis.shape[1:], vis.dtype)
		for i in xrange(len(self.nsplit)-1) : 
			if (self.verbose) : progressbar.Progress()
			self.gainnu[i] = vis[self.nsplit[i]:self.nsplit[i+1]].mean(0)
		#--------------------------------------------------
		verbose = self.masking.verbose
		self.masking.verbose = False
		self.gainnu, maskvalue = self.masking.MaskLoop(1, 7, 1, 3, 2, None, self.gainnu, None)
		self.masking.verbose = verbose
		self.gainnu.data[self.gainnu.mask] = maskvalue
		self.gainnu = jp.Smooth(self.gainnu.data, 1, 5, 1, Nprocess=self.Nprocess)
		# Normalize
		# The ratio of different gainnu is dependent on two things: (1) gain(t), (2) amplitudes of different regions of sky
	#	self.gainnu /= self.gainnu.mean(1)[:,None,:]
		#--------------------------------------------------
		# Plot
		titlehead = r'$g(\nu)$'
		fighead = 'gnu'
		x = self.antarray.Ant.freq
		xlabel = r'$\nu$ [MHz]'
		xaxes = plt_axes('x', 'both', [25,5])
		notNlinelist = None
		self._Plot(self.gainnu, 1, titlehead, fighead, outdir, x, xlabel, xaxes, Nline, notNlinelist, legendsize, legendloc, progressbar)
		if (self.verbose): print 'CaliGain.Gainnu:  end  @', jp.Time(1)+'\n'



#	def Gaint(self, timeper=3,timetimes=100, nfreq=None,nbl=None):
	def Gaint( self, timeper=3, timetimes=100, Nline=20, legendsize=10, legendloc=1, masktime=None ):
		'''
		Actually the return is:
			sigma_noise * gain(t)

		masktime:
			pair (n1,n2) or list of pair [(n1,n2), (n3,n4), ...]
			time range of bright source and so on
			For PAON4-CygA1dec, masktime=(5000,11000)

		nfreq, nbl:
			Plot a figure with nfreq, nbl

		nfreq:
			(1) int: freq bin
			(2) str with 'Hz': '1420.406MHz' or '1.420406GHz' or '1420406000Hz', etc
			(3) None: center bin

		nbl:
			(1) int: order of bl
			(2) int pair (int1, int2): channel pair of baseline
			(3) None: longest East-West baseline
		'''
		outdir = self.outdir + 'Gaint/'
		if (not os.path.exists(outdir)) : os.makedirs(outdir)
		#--------------------------------------------------
		if (masktime is not None) : 
			masktime = npfmt(masktime)
			if (masktime.size == 2) : masktime = masktime[None,:]
			masknsplit = np.zeros(len(self.nsplit)-1, bool)
			for i in xrange(len(self.nsplit)-1) : 
				n1, n2 = self.nsplit[i], self.nsplit[i+1]
				for j in xrange(len(masktime)) : 
					m1, m2 = np.sort(masktime[j])
					if(m1<=n1<=m2 or m1<=n2<=m2 or (n1<m1 and n2>m2)):
						masknsplit[i] = True
		#--------------------------------------------------
		if (self.verbose) : 
			print 'CaliGain.Gaint: start @', jp.Time(1)
			progressbar = jp.ProgressBar('    completed:', len(self.nsplit)-1+self.antarray.visorder.size)
		# Read data
		vistype = self.antarray.vistype[:-1]
		vis = self.antarray.vis[:,:,self.antarray.visorder]
		if (vistype == 'auto') : vis = vis.real
		if (len(vis.shape) == 2) : vis = vis[:,:,None]  # 3D
		# Replace masks
		vis[self.masking.mask] = self.masking.maskvalue
		# Get the fluctuation
		vis -= jp.Smooth(vis, 0, timeper, timetimes, Nprocess=self.Nprocess)
		# Masked
#		vis = np.ma.MaskedArray(vis, self.masking.mask)
		# Calculate the std of noise
		self.gaint = np.zeros((len(self.nsplit)-1,)+vis.shape[1:], vis.dtype)
		for i in xrange(len(self.nsplit)-1) : 
			if (self.verbose) : progressbar.Progress()
			if (vistype == 'auto') : self.gaint[i] = (vis[self.nsplit[i]:self.nsplit[i+1]]**2).mean(0)**0.5
			else : 
				self.gaint[i].real = (vis.real[self.nsplit[i]:self.nsplit[i+1]]**2).mean(0)**0.5
				self.gaint[i].imag = (vis.imag[self.nsplit[i]:self.nsplit[i+1]]**2).mean(0)**0.5
		#--------------------------------------------------
		# Reset the bad points
		verbose = self.masking.verbose
		self.masking.verbose = False
		self.gaint, maskvalue = self.masking.MaskLoop(0, 7, 1, 3, 2, None, self.gaint, None)
		self.masking.verbose = verbose
		mask = self.gaint.mask.copy()
		self.gaint.data[mask] = maskvalue
		self.gaint = jp.Smooth(self.gaint.data, 0, 5, 1, Nprocess=self.Nprocess)
		if (masktime is not None) : 
			self.gaint = np.ma.MaskedArray(self.gaint, np.zeros(self.gaint.shape, bool))
			self.gaint = self.gaint.T
			self.gaint.mask += masknsplit
			self.gaint = self.gaint.T
			mask = self.gaint.mask.copy()      #@
			maskvalue = self.gaint.data[mask]  #@
			self.gaint = jp.ResetMasked(self.gaint, 0, self.Nprocess)
		#--------------------------------------------------
		self.gaintmean = self.gaint / self.gaint.mean(0)
		n = jp.SelectLeastsq(self.gaintmean, 1, 10)
		self.gaintmean = self.gaint[:,n].mean(1)[:,None]
		if (masktime is not None) : self.gaint[mask] = maskvalue
		verbose = self.masking.verbose
		self.masking.verbose = False
		self.gaintmean, maskvalue = self.masking.MaskLoop(0, self.gaintmean.shape[0]/10, 5, 2, 2, None, self.gaintmean, None)
		self.masking.verbose = verbose
		self.gaintmean.data[self.gaintmean.mask] = maskvalue
		self.gaintmean = jp.Smooth(self.gaintmean.data, 0, self.gaintmean.shape[0]/20, 5, Nprocess=self.Nprocess)
		self.gaintmean = jp.Smooth(self.gaintmean, 0, 3, self.gaintmean.shape[0]/5, Nprocess=self.Nprocess)
		#--------------------------------------------------
		# Plot
		titlehead = r'$g(t)$'
		fighead = 'gaint'
		x = np.arange(self.gaint.shape[0])
		xlabel = 'time-axis'
		xaxes = None
		Nedge = self.gaint.shape[1]/10
		notNlinelist = range(0,Nedge)+range(self.gaint.shape[1]-Nedge,self.gaint.shape[1])
		self._Plot([self.gaint, self.gaintmean], 0, titlehead, fighead, outdir, x, xlabel, xaxes, Nline, notNlinelist, legendsize, legendloc, progressbar)
		if (self.verbose): print 'CaliGain.Gaint:  end  @', jp.Time(1)+'\n'



#
#
#
#
#
#
#
#
#
#
#
#		#---------- Plot ----------
#		nbl, strbl, norder = self._nbl(nbl)
#		nfreq, strfreq = self._nfreq(nfreq)
#		gaint = self.gaint[:,nfreq,norder]
#		# ylim
#		gaintm = spsn.medfilt(gaint.real, 7)
#		ymax, ymin = gaintm.max()*1.1, gaintm.min()*0.9
#		self._plot = (nbl,strbl,norder, nfreq,strfreq, ymax,ymin)
#		# Plot
#		x = np.arange(gaint.size)
#		if (vistype == 'auto') : 
#			plt.plot(x, gaint, 'bo', markersize=5, label='auto')
#			plt.legend()
#			plt.xlim(x.min(), x.max())
#			plt.ylim(ymin, ymax)
#			plt.xlabel('time', size=16)
#			plt.ylabel('[A.U.]', size=16)
#		else : 
#			plt.plot(x, gaint.real, 'bo', markersize=5,label='real')
#			plt.plot(x, gaint.imag, 'ro', markersize=5,label='imag')
#			plt.legend()
#			plt.xlim(x.min(), x.max())
#			plt.ylim(ymin, ymax)
#			plt.xlabel('time', size=16)
#			plt.ylabel('[A.U.]', size=16)
#		plt.title(r'$G(t) \cdot \sigma_n$, '+strfreq+', '+strbl, size=16)
#		plt.savefig(self.outdir+'CaliGain.Gaint_'+strfreq+'_'+strbl.split('=')[0]+'.png')
#		if (self.verbose) : print 'CaliGain.Gaint:  end  @', jp.Time(1)+'\n'
#
#
#
#	def Smooth( self, filtersize, smoothtimes ) : 
#		'''
#		Smooth self.gaint
#
#		filtersize:
#			int: use this int for both time-axis and freq-axis
#			int pair: (int1, int2), int1 for time-axis, int2 for freq-axis
#
#		smoothtimes:
#			int: for jp.Smooth() over time-axis
#		'''
#		if (self.verbose) : print 'CaliGain.Smooth: start @', jp.Time(1)
#		# filtersize
#		try : filtersize = filtersize[:2]
#		except : filtersize = (5, 3)
#		if (len(filtersize) == 1) : filtersize = (filtersize[0], filtersize[0])
#		# smoothtimes
#		try : smoothtimes = int(round(smoothtimes))
#		except : smoothtimes = self.gaint.shape[0]/5
#		# Do smoothing
#		vistype = self.antarray.vistype[:-1]
#		self.gaintsmooth = self.gaint*0
#		for i in xrange(self.gaint.shape[2]) : 
#			if (vistype == 'auto') : self.gaintsmooth[:,:,i] = spsn.medfilt2d(self.gaint[:,:,i], filtersize)
#			else : 
#				self.gaintsmooth[:,:,i].real = spsn.medfilt2d(self.gaint[:,:,i].real, filtersize)
#				self.gaintsmooth[:,:,i].imag = spsn.medfilt2d(self.gaint[:,:,i].imag, filtersize)
#		self.gaintsmooth = jp.Smooth(self.gaint, 0, 3, smoothtimes)
#		#--------------------------------------------------
#		#---------- Plot ----------
#		if ('_plot' in self.__dict__.keys()) : 
#			nbl,strbl,norder, nfreq,strfreq, ymax,ymin = self._plot
#			gaint  = self.gaint[:,nfreq,norder]
#			gaints = self.gaintsmooth[:,nfreq,norder]
#			# Plot
#			x = np.arange(gaint.size)
#			if (vistype == 'auto') : 
#				plt.plot(x, gaint, 'bo', markersize=3, label='auto')
#				plt.plot(x, gaints, 'm-', lw=3, label='auto smoothed')
#				plt.legend()
#				plt.xlim(x.min(), x.max())
#				plt.ylim(ymin, ymax)
#				plt.xlabel('time', size=16)
#				plt.ylabel('[A.U.]', size=16)
#			else : 
#				plt.plot(x, gaint.real, 'bo', markersize=3, label='real')
#				plt.plot(x, gaints.real, 'm-', lw=3, label='real smoothed')
#				plt.plot(x, gaint.imag, 'ro', markersize=3, label='imag')
#				plt.plot(x, gaints.imag, 'g-', lw=3, label='imag smoothed')
#				plt.legend()
#				plt.xlim(x.min(), x.max())
#				plt.ylim(ymin, ymax)
#				plt.xlabel('time', size=16)
#				plt.ylabel('[A.U.]', size=16)
#			plt.title(r'$G(t) \cdot \sigma_n$, '+strfreq+', '+strbl, size=16)
#			plt.savefig(self.outdir+'CaliGain.Gaint_'+strfreq+'_'+strbl.split('=')[0]+'.png')
#			if (self.verbose) : print 'CaliGain.Smooth:  end  @', jp.Time(1)+'\n'
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#	def See( self, timeper, timetimes, nfreq=None, nbl=None ) : 
#		'''
#		nfreq:
#			(1) int: freq bin
#			(2) str with 'Hz': '1420.406MHz' or '1.420406GHz' or '1420406000Hz', etc
#			(3) None: center bin
#
#		nbl:
#			(1) int: order of bl
#			(2) int pair (int1, int2): channel pair of baseline
#			(3) None: longest East-West baseline
#		'''
#		if (self.verbose) : print 'CaliGain.See: start @', jp.Time(1)
#		vistype = self.antarray.vistype[:-1]
#		nbl, strbl, norder = self._nbl(nbl)
#		nfreq, strfreq = self._nfreq(nfreq)
#		# time
#		timeper, timetimes = np.array([timeper, timetimes]).round().astype(int)
#		# Read data
#		vis = self.antarray.vis[:,:,self.antarray.visorder]
#		if (len(vis.shape) == 2) : vis = vis[:,:,None]
#		if (vistype == 'auto') : vis = vis.real
#		# Replace masks
#		vis[self.masking.mask] = self.masking.maskvalue
#		# Select
#		vis = vis[:,nfreq,norder]
#		# Get the fluctuation
#		viss = jp.Smooth(vis, 0, timeper, timetimes)
#		# Plot
#		x = np.arange(vis.size)
#		if (vistype == 'auto') : 
#			plt.figure(figsize=(8,8))
#			plt.subplot(2,1,1)
#			plt.plot(x, vis, 'b-', label='data')
#			plt.plot(x, viss, 'r-', label='smoothed')
#			plt.legend()
#			plt.xlabel('time', size=16)
#			plt.xlim(x.min(), x.max())
#			plt.ylabel('[A.U.]', size=16)
#			plt.subplot(2,1,2)
#			plt.plot(x, vis-viss, 'g-', label='residual')
#			plt.legend()
#			plt.xlabel('time', size=16)
#			plt.xlim(x.min(), x.max())
#			plt.ylabel('[A.U.]', size=16)
#		else : 
#			plt.figure(figsize=(17,8))
#			plt.subplot(2,2,1)
#			plt.plot(x, vis.real, 'b-', label='data.real')
#			plt.plot(x, viss.real, 'r-', label='smoothed.real')
#			plt.legend()
#			plt.xlabel('time', size=16)
#			plt.xlim(x.min(), x.max())
#			plt.ylabel('[A.U.]', size=16)
#			plt.subplot(2,2,3)
#			plt.plot(x, vis.real-viss.real, 'g-', label='residual.real')
#			plt.legend()
#			plt.xlabel('time', size=16)
#			plt.xlim(x.min(), x.max())
#			plt.ylabel('[A.U.]', size=16)
#			plt.subplot(2,2,2)
#			plt.plot(x, vis.imag, 'b-', label='data.imag')
#			plt.plot(x, viss.imag, 'r-', label='smoothed.imag')
#			plt.legend()
#			plt.xlabel('time', size=16)
#			plt.xlim(x.min(), x.max())
#			plt.ylabel('[A.U.]', size=16)
#			plt.subplot(2,2,4)
#			plt.plot(x, vis.imag-viss.imag, 'g-', label='residual.imag')
#			plt.legend()
#			plt.xlabel('time', size=16)
#			plt.xlim(x.min(), x.max())
#			plt.ylabel('[A.U.]', size=16)
#		plt.suptitle('Flucation, '+strfreq+', '+strbl, size=20)
#		plt.savefig(self.outdir+'CaliGain.See_'+strfreq+'_'+strbl.split('=')[0]+'.png')
#		plt.close()
#		if (self.verbose): print 'CaliGain.See:  end  @', jp.Time(1)+'\n'
#

	def _Plot( self, arraypair, axis, titlehead, fighead, outdir, xlist=None, xlabel='', xaxes=None, Nline=20, notNlinelist=None, legendsize=10, legendloc=1, progressbar=None ) : 
		'''
		arraypair: 
			(1) list or tuple: array, arraymean = arraypair
			(2) ndarray: array, arraymean = arraypair, None

		axis:
			Which axis to be plot
			If want to plot freq-Amp, axis=1
			Else want to plot time-Amp, axis=0

		titlehead:
			plt.title(r''+ titlehead +' of channel '+str(chanidx)+', total '+str(len(color))+' curves')

		fighead:
			plt.savefig( fighead +'_chan'+str(chanidx)+'_N'+str(len(color))+'.png')

		xlist:
			If =None, set xlist=np.arange(N)

		xlabel:
			plt.xlabel(xlabel)

		xaxes:
			xaxes = plt_axes('x', 'both', [25,5])

		Nline:
			How many lines do you want to plot

		notNlinelist:
			In this function, we will get Nlinelist, if Nlinelist[i] in notNlinelist, throw it (don't plot)

		legendsize:
			Use for plt.legend(fontsize=legendsize)
		'''
		istype = jp.IsType()
		if (istype.islist(arraypair) or istype.istuple(arraypair)):
			array, arraymean = arraypair
			arraymean /= arraymean.mean(0)
		else : array, arraymean = arraypair, None
		try : Nline = int(round(Nline))
		except : Nline = 20
		axis = int(round(axis))
		if (axis not in [0,1]) : jp.Raise('axis must be 0=>time OR 1=>freq')
		axis2 = 0 if(axis==1) else 1
		if(Nline>array.shape[axis2]):Nline=self.gainnu.shape[axis2]
		try : notNlinelist = list(notNlinelist)
		except : notNlinelist = []
		#--------------------------------------------------
		arraysmooth =jp.Smooth(array,axis,array.shape[axis]/10, 1)
		darray = (array - arraysmooth).std(axis)
		#--------------------------------------------------
		for i in xrange(array.shape[2]) : 
			if (self.verbose) : progressbar.Progress()
			chanidx = self.antarray.Blorder.blorder[self.antarray.visorder][i][0]
			#------------------- Get Nline --------------------
			nsplit = np.linspace(0, array.shape[axis2], Nline+1).astype(int)
			darrayi = darray[:,i]+1j*np.arange(darray[:,i].size) 
			Nlinelist = []
			for k in xrange(Nline) : 
				darrayk = darrayi[nsplit[k]:nsplit[k+1]]
				darrayk = np.sort(darrayk).imag.astype(int)[0]
				if (darrayk not in notNlinelist) : Nlinelist.append(darrayk)
			Nlinelist, darrayi = np.array(Nlinelist), 0
			color = plt_color(len(Nlinelist))
			#--------------------------------------------------
			if (axis == 1) : arrayi = array[Nlinelist,:,i]
			else : arrayi = array[:,Nlinelist,i]
			arrayimean = arrayi.mean(axis)
			if (axis == 1) : arrayimean = arrayimean[:,None]
			arrayi /= arrayimean
			vmin, vmax = arrayi.min(axis), arrayi.max(axis)
			vmin, vmax = vmin.max(), vmax.min()
			vmin, vmax = vmin-0.3*(vmax-vmin), vmax+0.3*(vmax-vmin)
			if (vmin < 0) : vmin = 0
			for j in xrange(len(Nlinelist)) : 
				if (axis == 1) : plt.plot(xlist, arrayi[j], color=color[j], label=str(Nlinelist[j]+1))
				else : plt.plot(xlist, arrayi[:,j], color=color[j], label=str(Nlinelist[j]+1))
			if (arraymean is not None) : 
				plt.plot(xlist, arraymean[:,i], 'k-', lw=3, label='Guess')
			plt.legend(fontsize=legendsize, loc=legendloc)
			plt.xlabel(xlabel, size=16)
			plt.xlim(int(round(xlist.min())), int(round(xlist.max())))
			xaxes
			plt.ylabel('[A.U.]', size=16)
			plt.ylim(vmin, vmax)
			if (arraymean is None) : 
				plt.title(titlehead+' of channel '+str(chanidx)+', total '+str(len(color))+' curves', size=16)
				plt.savefig(outdir+fighead+'_chan'+str(chanidx)+'_N'+str(len(color))+'.png')
			else : 
				plt.title(titlehead+' of channel '+str(chanidx)+', total '+str(len(color))+'+1 curves', size=16)
				plt.savefig(outdir+fighead+'_chan'+str(chanidx)+'_N'+str(len(color))+'-1.png')
			plt.close()



#
#
#
#	def _nfreq( self, nfreq=None ) : 
#		'''
#		nfreq:
#			(1) int: freq bin
#			(2) str with 'Hz': '1420.406MHz' or '1.420406GHz' or '1420406000Hz', etc
#			(3) None: center bin
#
#		return [nfreq, strfreq]
#			nfreq: AntArray.Ant.freq[nfreq], AntArray.vis[:,nfreq]
#		'''
#		istype = jp.IsType()
#		try : 
#			nfreq = int(round(nfreq))
#			strfreq = '%.1fMHz' % self.antarray.Ant.freq[nfreq]
#		except : 
#			if (istype.isstr(nfreq) and 'hz' in nfreq.lower()) : 
#				strfreq = nfreq
#				n = nfreq.lower().rfind('hz')
#				if (nfreq[n-1].lower() == 'm') : nfreq = float(nfreq[:n-1])  # MHz
#				elif (nfreq[n-1].lower() == 'g') : nfreq = float(nfreq[:n-1])*1000
#				else : nfreq = float(nfreq[:n-1])*1e-6
#				nfreq = abs(self.antarray.Ant.freq-nfreq)
#				nfreq = np.where(nfreq==nfreq.min())[0][0]
#			else : 
#				nfreq = self.antarray.vis.shape[1] /2
#				strfreq = '%.1fMHz' % self.antarray.Ant.freq[nfreq]
#		return [nfreq, strfreq]
#
#
#	def _nbl( self, nbl=None ) : 
#		'''
#		nbl:
#			(1) int: order of bl
#			(2) int pair (int1, int2): channel pair of baseline
#			(3) None: longest East-West baseline
#		return [nbl, strbl, norder]
#
#		nbl: AntArray.Blorder.blorder[nbl], AntArray.vis[:,:,nbl]
#		norder: AntArray.visorder[norder], AntArray.vis[:,:,AntArray.visorder[norder]]
#		'''
#		vistype = self.antarray.vistype[:-1]
#		try : 
#			nbl = int(round(nbl))
#			if (nbl not in self.antarray.visorder) : raise
#		except : 
#			try : 
#				nbl = self.antarray.Blorder.Bl2Order(nbl[:2])
#				if (nbl not in self.antarray.visorder) : raise
#			except : 
#				bl = abs(self.antarray.Blorder.baseline[self.antarray.visorder][:,0])
#				bl = np.where(bl==bl.max())[0][0]
#				nbl = self.antarray.visorder[bl]
#		norder = self.antarray.Blorder.blorder[nbl]
#		bl     = self.antarray.Blorder.baseline[nbl]
#		if (vistype == 'auto') : strbl = 'auto:%i' % norder[0]
#		else : strbl = 'cross:%i-%i=(%.3f, %.3f, %.3f)' % (tuple(norder)+tuple(bl))
#		norder = np.where(self.antarray.visorder==nbl)[0][0]
#		return [nbl, strbl, norder]
#
#
