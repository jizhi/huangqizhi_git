import os
import numpy as np
import jizhipy as jp
from jizhipy.Plot import *
import scipy.signal as spsn
##################################################



class CaliGain( object ) : 


	def __init__( self, antarray=None, masking=None, Nprocess=None, verbose=True, outdir=None ) : 
		class _Params( object ) : pass
		self.Params = _Params()
		self.Params.init = {'Nprocess':Nprocess, 'verbose':verbose}
		self.Nprocess, self.verbose = Nprocess, verbose
		if (antarray is not None) : self.antarray = antarray
		if (masking  is not None) : self.masking  = masking
		self.outdir = jp.Outdir((None,'file'), (0,'file'))
		if (outdir is not None) : self.outdir = jp.Mkdir(self.outdir+outdir)
		self.starttime = jp.Time(1)


	def Window( self, size ) : 
		''' How many time points to calculate one std? '''
		self.Params.Window = {'size':size}
		N = int(round(self.antarray.vis.shape[0] / size))
		self.nsplit = np.linspace(0, self.antarray.vis.shape[0], N+1).round().astype(int)



	def Gainnu( self, Nline=20, gainnuper=3, gainnutimes=10, legendsize=9, legendloc=1 ) : 
		''' Calculate g(nu) and plot '''
		self.Params.Gainnu = {'Nline':Nline, 'legendsize':legendsize, 'legendloc':legendloc}
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliGain.Gainnu: start @', starttime
			progressbar = jp.ProgressBar('    completed:', len(self.nsplit)-1+self.antarray.visorder.size+1)
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
		self.masking.verbose = False
		self.gainnu, maskvalue = self.masking.MaskLoop(1, 7, 1, 3, 2, None, self.gainnu, None)
		self.masking.verbose = self.masking.Params.init['verbose']
		self.gainnu.data[self.gainnu.mask] = maskvalue
		self.gainnu = jp.Smooth(self.gainnu.data, 1, 5, 1, Nprocess=self.Nprocess)
		# Normalize
		# The ratio of different gainnu is dependent on two things: (1) gain(t), (2) amplitudes of different regions of sky
	#	self.gainnu /= self.gainnu.mean(1)[:,None,:]
		# Normalize
		gainnu = self.gainnu / self.gainnu.mean(1)[:,None,:]
		self.gainnumean, Nlinelist = self._Gainmean(gainnu, 1, Nline, None, gainnuper, gainnutimes, False)
		np.save(self.outdir+'gainnu.npy', self.gainnu)
		np.save(self.outdir+'Nlinelist_nu.npy', Nlinelist)
		np.save(self.outdir+'gainnumean.npy', self.gainnumean)
		#--------------------------------------------------
		# Plot
		titlehead = r'$g(\nu)$'
		fighead = 'gnu'
		x = self.antarray.Ant.freq
		xlabel = r'$\nu$ [MHz]'
		xaxes = plt_axes('x', 'both', [25,5])
		self._Plot(gainnu[Nlinelist], self.gainnumean, 1, titlehead, fighead, self.outdir, x, xlabel, xaxes, Nlinelist, legendsize, legendloc, progressbar)
		if (self.verbose) : 
			print '    gainnu.npy, gainnumean.npy, Nlinelist.npy  saved'
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliGain.Gainnu:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'



	def Gaint( self, timeper=3, timetimes=100, gaintper=3, gainttimes=100, Nline=20, legendsize=10, legendloc=1, masktime=None ) :
		'''
		timeper, timetimes:
			vis - Smooth(vis, 0, timeper, timetimes)

		gaintper, gainttimes:
			Smooth(gaintmean, 0, gaintper, gainttimes)

		Actually the return is:
			sigma_noise * gain(t)

		Nline: 
			How many curves used to calculate gaintmean and be plotted

		legendsize, legendloc:
			plt.legend(fontsize=legendsize, loc=legendloc)

		masktime:
			pair (n1,n2) or list of pair [(n1,n2), (n3,n4), ...]
			time range of bright source and so on
			For PAON4-CygA1dec, masktime=(5000,11000)
		'''
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliGain.Gaint: start @', starttime
			progressbar = jp.ProgressBar('    completed:', len(self.nsplit)-1+self.antarray.visorder.size+5)
		#--------------------------------------------------
		self.Params.Gaint = {'timeper':timeper, 'timetimes':timetimes, 'gaintper':gaintper, 'gainttimes':gainttimes, 'Nline':Nline, 'legendsize':legendsize, 'legendloc':legendloc, 'masktime':masktime}
		#--------------------------------------------------
		# Read data
		vistype = self.antarray.vistype[:-1]
		vis = self.antarray.vis[:,:,self.antarray.visorder]
		if (vistype == 'auto') : vis = vis.real
		if (len(vis.shape) == 2) : vis = vis[:,:,None]  # 3D
		if (self.verbose) : progressbar.Progress()
		#--------------------------------------------------
		# Replace masks
		vis[self.masking.mask] = self.masking.maskvalue
		# Get the fluctuation
		vis -= jp.Smooth(vis, 0, timeper, timetimes, Nprocess=self.Nprocess)
		# Masked
#		vis = np.ma.MaskedArray(vis, self.masking.mask)
		if (self.verbose) : progressbar.Progress()
		#--------------------------------------------------
		# Calculate the std of noise
		# self.gaint: auto=>real, cross=>complex
		self.gaint = np.zeros((len(self.nsplit)-1,)+vis.shape[1:], vis.dtype)
		for i in xrange(len(self.nsplit)-1) : 
			if (self.verbose) : progressbar.Progress()
			if (vistype == 'auto') : self.gaint[i] = (vis[self.nsplit[i]:self.nsplit[i+1]]**2).mean(0)**0.5
			else : 
				self.gaint[i].real = (vis.real[self.nsplit[i]:self.nsplit[i+1]]**2).mean(0)**0.5
				self.gaint[i].imag = (vis.imag[self.nsplit[i]:self.nsplit[i+1]]**2).mean(0)**0.5
		#--------------------------------------------------
		#--------------------------------------------------
		# Mask the bad/RFI points
		if (self.verbose) : progressbar.Progress()
		self.masking.verbose = False
		self.gaint, maskvalue = self.masking.MaskLoop(0, self.gaint.shape[0]/10, 1, 3, 2, None, self.gaint, None)  #@#@
		self.masking.verbose = self.masking.Params.init['verbose']
		self.gaint.data[self.gaint.mask] = maskvalue
		self.gaint = jp.Smooth(self.gaint.data, 0, 5, 1, Nprocess=self.Nprocess)  # Use this as the original gaint, not Masked
		#--------------------------------------------------
		# Mask manually if masktime
		# masktime => masknsplit (bool)
		if (self.verbose) : progressbar.Progress()
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
			#----------
			self.gaint = np.ma.MaskedArray(self.gaint, np.zeros(self.gaint.shape, bool))
			self.gaint = self.gaint.T
			self.gaint.mask += masknsplit
			self.gaint = self.gaint.T
			gaintreset = jp.ResetMasked(self.gaint, 0, self.Nprocess)
			self.gaint = self.gaint.data
		else : gaintreset = self.gaint
		# self.gaint, gaintreset
		# Normalize self.gaint
		gaint = self.gaint.copy()
		gm = gaintreset.mean(0)
		if (vistype == 'auto') : 
			gaint /= gaint.mean(0)
			gaintreset /= gm
		else : 
			gaint.real /= gaint.real.mean(0)
			gaint.imag /= gaint.imag.mean(0)
			gaintreset.real /= gm.real
			gaintreset.imag /= gm.imag
		#--------------------------------------------------
		Nedge = self.gaint.shape[1]/10
		notNlinelist = range(0,Nedge)+range(self.gaint.shape[1]-Nedge,self.gaint.shape[1])
		self.gaintmean, Nlinelist = self._Gainmean(gaintreset, 0, Nline, notNlinelist, gaintper, gainttimes, True)
		#--------------------------------------------------
		# self.norm = self.gaint.mean(0) / self.gaintmean.mean(0)
		gmm = self.gaintmean.mean(0)
		if (vistype == 'auto') : self.norm = (gm / gmm)[None,:]
		else : 
			normr = (gm.real / gmm.real)[None,:]
			normi = (gm.imag / gmm.imag)[None,:]
			self.norm = normr + 1j*normi
		gm = gmm = 0 #@
		#--------------------------------------------------
		np.save(self.outdir+'gaint.npy', self.gaint)
		np.save(self.outdir+'Nlinelist_t.npy', Nlinelist)
		np.save(self.outdir+'gaintmean.npy', self.gaintmean)
		# self.gaintmean is normalized and self.gaint is NOT !
		#--------------------------------------------------
		# Plot
		titlehead = r'$g(t)$'
		fighead = 'gaint'
		x = np.arange(self.gaint.shape[0])
		xlabel = 'time-axis'
		xaxes = None
		self._Plot(gaint[:,Nlinelist], self.gaintmean, 0, titlehead, fighead, self.outdir, x, xlabel, xaxes, Nlinelist, legendsize, legendloc, progressbar)
		if (self.verbose) : 
			print '    gaint.npy, gaintmean.npy, Nlinelist.npy  saved'
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliGain.Gaint:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'



	def _Nlinelist( self, array, axis, Nline, notNlinelist=[] ) : 
		'''
		Nline:
			How many lines do you want to plot
			For example, plot Gaint, there are 512 freq points, total 512 lines. Set Nline=20, will just plot 20 (uniform interval) of them
			Will get Nlinelist (list of int/index)

		notNlinelist:
			Parameter Nline will get Nlinelist, if Nlinelist[i] in notNlinelist, throw it (don't plot)
		'''
		axis = int(round(axis))
		if (axis not in [0,1]) : jp.Raise('axis must be 0=>time OR 1=>freq')
		# Move axis to 0
		array = jp.ArrayAxis(array, axis, 0, 'move')
		#--------------------------------------------------
		# Number of curves in the figure
		try : Nline = int(round(Nline))
		except : Nline = 20
		# time or freq to be plotted
		if(Nline > array.shape[1]) : Nline = array.shape[1]
		#--------------------------------------------------
		try : notNlinelist = list(notNlinelist)
		except : notNlinelist = []
		#--------------------------------------------------
		# Judge which slice to be plotted
		darray = (array - jp.Smooth(array, 0, array.shape[0]/10, 1)).std(0)  # real value, 2D, n1xn2
		# Select leastsq baseline
		ddarray = darray.std(0)
		ddarray = ddarray + 1j*np.arange(ddarray.size)
		nbl = np.sort(ddarray).imag.astype(int)[0]
		darray = darray[:,nbl]
		ddarray = 0 #@
		#--------------------------------------------------
		Nlinelist = []
		nsplit= np.linspace(0, array.shape[1], Nline+1).astype(int)
		darray = darray + 1j*np.arange(darray.size) 
		for i in xrange(Nline) : 
			n = darray[nsplit[i]:nsplit[i+1]]
			n = np.sort(n).imag.astype(int)[0]
			if (n not in notNlinelist) : Nlinelist.append(n)
		return array, axis, np.array(Nlinelist)



	def _Gainmean( self, array, axis, Nline, notNlinelist, gainper, gaintimes, mask=False ) : 
		array, axis, Nlinelist = self._Nlinelist(array, axis, Nline, notNlinelist)
		array = array[:,Nlinelist,:]
		if (mask) : 
			self.masking.verbose = False
			array, maskvalue = self.masking.MaskLoop(0, array.shape[0]/10, 1, 3, 1, None, array, None)
			self.masking.verbose=self.masking.Params.init['verbose']
			array.data[array.mask] = maskvalue
			array = array.data.mean(1)
		else : array = array.mean(1)
		array = jp.Smooth(array, 0, gainper, gaintimes)
		# Normalize again
		if (array.dtype.name[:3] != 'com') : array /= array.mean(0)
		else : 
			array.real /= array.real.mean(0)
			array.imag /= array.imag.mean(0)
		if (axis == 0) : array = array[:,None,:]
		else : array = array[None,:,:]
		return array, Nlinelist
		


	def _Plot( self, array, arraymean, axis, titlehead, fighead, outdir, xlist=None, xlabel='', xaxes=None, Nlinelist=None, legendsize=10, legendloc=1, progressbar=None ) : 
		'''
		If plot, need to normalize, because different freq gets different scale of gaint

		array, arraymean:
			array to be plotted

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

		legendsize:
			Use for plt.legend(fontsize=legendsize)
		'''
		if (self.verbose) : progressbar.Progress()
		if (arraymean.dtype.name[:3] != 'com') : vistype = 'auto'
		else : vistype = 'cross'
		axis = int(round(axis))
		if (axis not in [0,1]) : jp.Raise('axis must be 0=>time OR 1=>freq')
		array = jp.ArrayAxis(array, axis, 0, 'move')
		arraymean = jp.ArrayAxis(arraymean, axis, 0, 'move')
		#--------------------------------------------------
		for i in xrange(arraymean.shape[2]) : # Plot each baseline
			ai, ami = array[:,:,i], arraymean[:,0,i]
			if (self.verbose) : progressbar.Progress()
			chanidx = self.antarray.Blorder.blorder[self.antarray.visorder][i]
			if (chanidx[0]==chanidx[1]): chanidxstr= str(chanidx[0])
			else : chanidxstr = str(chanidx[0])+'-'+str(chanidx[1])
			#--------------------------------------------------
			if (vistype=='auto') : color = plt_color(array.shape[1])
			else : color = plt_color(2*array.shape[1])
			#--------------------------------------------------
			# vmin, vmax
			def _MinMax( array, ratio=0.5 ) : 
				vmin, vmax = array.min(), array.max()
				vmin, vmax = vmin-ratio*(vmax-vmin), vmax+ratio*(vmax-vmin)
				if (vmin < 0) : vmin = 0
				return vmin, vmax
			ratio = 0.5 if(axis==0)else 0.2
			if (vistype == 'auto') : vmin, vmax = _MinMax(ami,ratio)
			else : 
				vminr, vmaxr = _MinMax(ami.real, ratio)
				vmini, vmaxi = _MinMax(ami.imag, ratio)
				vmin, vmax = min(vminr, vmini), max(vmaxr, vmaxi)
			#--------------------------------------------------
			if (axis == 1) : lw = np.linspace(5, 1, ai.shape[1]+1)
			for j in xrange(ai.shape[1]) : 
				# Plot array
				if (vistype == 'auto') :
					if (axis == 0) : plt.plot(xlist, ai[:,j], color=color[j], ls='', marker='o', markersize=3, label=str(Nlinelist[j]+1))
					else : plt.plot(xlist, ai[:,j], color=color[j], label=str(Nlinelist[j]+1), lw=lw[j])
				else : 
					if (axis == 0) : 
						plt.plot(xlist, ai[:,j].real, color=color[2*j+1], ls='', marker='o', markersize=3, label=str(Nlinelist[j]+1)+', real')
						plt.plot(xlist, ai[:,j].imag, color=color[2*j+1], ls='', marker='o', markersize=3, label=str(Nlinelist[j]+1)+', imag')
					else : 
						plt.plot(xlist, ai[:,j].real, color=color[2*j+1], label=str(Nlinelist[j]+1)+', real', lw=lw[j])
						plt.plot(xlist, ai[:,j].imag, color=color[2*j+1], label=str(Nlinelist[j]+1)+', imag', lw=lw[j])
			# Plot arraymean
			lw = 3 if(axis==0)else 1
			if (vistype == 'auto') : 
				plt.plot(xlist, ami, 'b-', lw=lw, label='mean')
			else : 
				plt.plot(xlist, ami.real, 'c-', lw=lw, label='mean, real')
				plt.plot(xlist, ami.imag, 'b-', lw=lw, label='mean, imag')
			plt.legend(fontsize=legendsize, loc=legendloc)
			plt.xlabel(xlabel, size=16)
			plt.xlim(int(round(xlist.min())),int(round(xlist.max())))
			xaxes
			plt.ylabel('[A.U.]', size=16)
			plt.ylim(vmin, vmax)
			nstr = '1' if(vistype=='auto')else '2'
			plt.title(titlehead+' of channel '+chanidxstr+', total '+str(ai.shape[1])+'+'+nstr+' curves', size=16)
			plt.savefig(outdir+fighead+'_chan'+chanidxstr+'_N'+str(ai.shape[1])+'.png')
			plt.close()


