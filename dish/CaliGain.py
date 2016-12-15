import jizhipy as jp
from jizhipy.Plot import *
import numpy as np

from Masking import *
##################################################





class CaliGain( object ) : 



	def __init__( self, masking=None, verbose=True, outdir='' ) : 
		try : 
			self.verbose = masking.verbose
			self.masking = masking
		except : 
			self.verbose = bool(verbose)
			self.masking = Masking(verbose=False)
			self.masking.verbose = self.verbose
		if (self.verbose) : print '-------------------- CaliGain --------------------\n'
		self.outdir = jp.Outdir((None,'file'), (0,'file'))
		try : self.outdir = jp.Mkdir(self.outdir + outdir)
		except : pass
		#--------------------------------------------------
		self.Nprocess = self.masking.Nprocess
		self.antarray = self.masking.antarray  # class:AntArray | np.ndarray
		self.antarray.verbose = self.masking.verbose
		try : self.timeper, self.timetimes = self.masking.per, self.masking.times
		except : self.timeper, self.timetimes = 3, 100
		try : self.outname = self.outdir + 'CaliGain_'+self.antarray.File.filedirname+'.hdf5'
		except : pass
		#--------------------------------------------------
		self.nsplit = None
		self.gainnuper, self.gainnutimes, self.gainnu, self.freq21cm, self.gainnu21cm, self.gainnumean, self.Nlinelist_gainnu = None, None, None, None, None, None, None
		self.gaintper, self.gainttimes, self.gaint, self.gaintmean, self.Nlinelist_gaint = None, None, None, None, None





	def Window( self, size ) : 
		''' How many time points to calculate one std? '''
		N = int(round(self.antarray.vis.shape[0] / size))
		self.nsplit = np.linspace(0, self.antarray.vis.shape[0], N+1).round().astype(int)
		if (self.verbose) : 
			print 'Masking.Window'
			print '    windowsize='+str(int(round(size)))
			print '    self.nsplit.size='+str(self.nsplit.size)+'\n'





	def Gainnu( self, gainnuper=3, gainnutimes=10, Nline=20, legendsize=9, legendloc=1 ) : 
		''' Calculate g(nu) and plot '''
		if (self.verbose) : 
			print 'CaliGain.Gainnu'
			print '    start @', jp.Time(1)
			progressbar = jp.ProgressBar('    completed:', len(self.nsplit)-1+self.antarray.vis.shape[2]+1)
		vistype = 'cross' if(self.antarray.vis.dtype.name[:7]=='complex')else 'auto'
		if (vistype != 'auto') : 
			print '    Warning: should use "auto-correlation", but now AntArray.vistype="'+vistype+'"'
			print '    Do nothing !'
			if (self.verbose): print '    end   @', jp.Time(1)+'\n'
			return
		#--------------------------------------------------
		if (self.antarray.vis.shape[1] == 1) : 
			print '    Warning: self.antarray.shape[1] == 1'
			print '    Can NOT do CaliGain.Gainnu() !'
			if (self.verbose): print '    end   @', jp.Time(1)+'\n'
			return
		#--------------------------------------------------

		self.gainnuper, self.gainnutimes = int(round(gainnuper)), int(round(gainnutimes))
		vis=np.ma.MaskedArray(self.antarray.vis, self.masking.mask)
		# Calculate
		self.gainnu = np.zeros((self.nsplit.size-1,)+vis.shape[1:], vis.data.dtype)
		for i in xrange(len(self.nsplit)-1) : 
			if (self.verbose) : progressbar.Progress()
			self.gainnu[i] = vis[self.nsplit[i]:self.nsplit[i+1]].mean(0)
		#--------------------------------------------------

		# Save around 1420.40575 MHz to see the Doppler shift
		try : 
			freq = self.antarray.freq
			tf = (1417.5<=freq)*(freq<=1423.5)
			self.freq21cm = freq[tf]
			if (self.freq21cm.size == 0) : print "    Warning: self.antarray.freq does NOT include 1417.5~1423.5MHz, NOT save fo['gainnu21cm'], fo['freq21cm'"
			else : self.gainnu21cm = self.gainnu[:,tf,:].copy()
		except : pass
		#--------------------------------------------------

		self.gainnu = jp.Medfilt(self.gainnu, 1, 3, self.Nprocess)
		# Remove bad value from gainnu
		self.gainnu = jp.ArrayAxis(self.gainnu, 1, 0, 'move')
		antarray = AntArray(verbose=False)
		antarray.vis = self.gainnu
		masking = Masking(antarray)
		masking.MaskLoop(0,self.gainnuper, self.gainnutimes, 4,5,0)
		self.gainnu = jp.ArrayAxis(masking.antarray.vis,0,1,'move')
		self.gainnu = jp.Smooth(self.gainnu, 1, 5, 1, Nprocess=self.Nprocess)
		#--------------------------------------------------

		# Normalize
		# The ratio of different gainnu is dependent on two things: (1) gain(t), (2) amplitudes of different regions of sky
	#	self.gainnu /= self.gainnu.mean(1)[:,None,:]
		self.gainnumean, self.Nlinelist_gainnu, arrstd = self._Gainmean(self.gainnu, 1, Nline, None, gainnuper, gainnutimes)
		#--------------------------------------------------

		# Plot
		titlehead = r'$g(\nu)$'
		fighead = 'gnu'
		xlabel = r'$\nu$ [MHz]'
		xaxes = ('x', 'both', [25,5])
		gainnu = self.gainnu / self.gainnu.mean(1)[:,None,:]
		gainnu /= self.gainnumean.mean(1)[:,None,:]
		self._Plot(gainnu[self.Nlinelist_gainnu], self.gainnumean, 1, titlehead, fighead, self.outdir, freq, xlabel, xaxes, self.Nlinelist_gainnu, legendsize, legendloc, progressbar)
		if (self.gainnu21cm is not None) : self._Plot21cm()
		if (self.verbose) : print '    end   @', jp.Time(1)+'\n'





	def _Plot21cm( self ) : 
		if (self.verbose) : '    Plotting  gnu_21cm.png'
		freq = self.freq21cm.copy()
		tf = (1419<=freq)*(freq<=1422)
	#	gnu21 = self.gainnu21cm[self.Nlinelist_gainnu,tf]
		gnu21 = self.gainnu21cm[self.Nlinelist_gainnu][:,tf]
		gnu21 /= (gnu21[:,:1] + gnu21[:,-1:])/2.
		gnu21 = gnu21.sum(1)
		gainnu21cm = []
		for i in xrange(gnu21.shape[-1]) : 
			n = np.where(gnu21[:,i]==gnu21[:,i].max())[0][0]
			gainnu21cm.append( self.gainnu21cm[self.Nlinelist_gainnu[n],:,i] )
		nv = np.arange(len(gainnu21cm))
		nf, strfreq, nv, strblo, strbli, blname = self.antarray.Plotnfnv(0, nv, True)
		blname = 'chan'
	#	s = nv.size**0.5
	#	if (s == int(s)) : s = int(s)
	#	else : s = int(s)+1
	#	if (s*(s-1) >= nv.size) : sr, sc = s-1, s
	#	else : sr = sc = s
		sr, sc = int(round(len(nv)/2.+1e-3)), 2
		plt.figure(figsize=(sc*16/3., sr*4))
		for i in xrange(len(nv)) : 
			plt.subplot(sr, sc, i+1)
			a = gainnu21cm[i]
			plt.plot(self.freq21cm, a, 'b-')
			plt.plot(self.freq21cm, a, 'ro', markersize=3, label=blname+'\n '+strblo[i])
			plt.legend(fontsize=12)
			plt.xlim(self.freq21cm[0], self.freq21cm[-1])
			plt_axes('x', 'both', [1, 0.1])
		#	if (i % 4 == sr-1) : plt.xlabel(r'$\nu$ [MHz]')
			plt.xlabel(r'$\nu$ [MHz]', size=16)
		#	if (i % sc == 0) : plt.ylabel('[A.U.]')
			d = a.max()-a.min()
			plt.ylim(a.min()-d*0.05, a.max()+d*0.1)
		plt.suptitle(r'$g(\nu)$ around 1420.40575MHz', size=16)
		plt.savefig(self.outdir+'gnu_21cm.png')
		plt.close()





	def Gaint( self, timeper, timetimes, gaintper=3, gainttimes=100, Nline=20, legendsize=10, legendloc=1, masktime=None):
		'''
		timeper, timetimes:
			vis - Smooth(vis, 0, timeper, timetimes) to calculate the sigma of noise

		gaintper, gainttimes:
			Smooth(gaintmean, 0, gaintper, gainttimes)

		Actually the return is:
			sigma_noise * gain(t)

		Nline: 
			How many curves used to calculate gaintmean and be plotted

		legendsize, legendloc:
			plt.legend(fontsize=legendsize, loc=legendloc)

		masktime:
			NOTE THAT: masktime should be obtained from Auto-corr, not from Cross-coor, because in Cross, can't see the diffuse sources, but these diffuse sources also affect gaint ! However, we can see these diffuse sources clearly in Auto.
			pair (n1,n2) or list of pair [(n1,n2), (n3,n4), ...]
			time range of bright source and so on
			For PAON4-CygA1dec, masktime=(5000,11000)
		'''
		if (self.verbose) : 
			print 'CaliGain.Gaint'
			print '    start @', jp.Time(1)
			progressbar = jp.ProgressBar('    completed:', len(self.nsplit)-1+self.antarray.visorder.size+1)
		if (timeper   is None) : timeper   = self.timeper
		if (timetimes is None) : timetimes = self.timetimes
		self.timeper, self.timetimes = int(round(timeper)), int(round(timetimes))
		self.gaintper, self.gainttimes = int(round(gaintper)), int(round(gainttimes))
		#--------------------------------------------------

		# Get the fluctuation
		if (timeper < 3) : timeper = 3
		vis = self.antarray.vis.copy()
		vis[self.masking.mask] = self.masking.maskvalue
		vis = vis - jp.Smooth(vis, 0, timeper, timetimes, Nprocess=self.Nprocess)
		# Masked
		vis = np.ma.MaskedArray(vis, self.masking.mask)
		#--------------------------------------------------

		# Calculate the std of noise
		# self.gaint: auto=>real, cross=>complex
		self.gaint = np.zeros((len(self.nsplit)-1,)+vis.shape[1:], vis.dtype)
		vistype = 'cross' if(vis.dtype.name[:7]=='complex')else 'auto'
		for i in xrange(len(self.nsplit)-1) : 
			if (self.verbose) : progressbar.Progress()
			if (vistype == 'auto') : self.gaint[i] = (vis[self.nsplit[i]:self.nsplit[i+1]]**2).mean(0)**0.5
			else : 
				self.gaint[i].real = (vis.real[self.nsplit[i]:self.nsplit[i+1]]**2).mean(0)**0.5
				self.gaint[i].imag = (vis.imag[self.nsplit[i]:self.nsplit[i+1]]**2).mean(0)**0.5
		#--------------------------------------------------
		#--------------------------------------------------

		self.gaint = jp.Medfilt(self.gaint, 0, 3)
		# Mask the bad/RFI points
		antarray = AntArray(verbose=False)
		antarray.vis = self.gaint
		masking = Masking(antarray)
		masking.MaskLoop(0, self.gaintper, self.gainttimes, 4,5,0)
		self.gaint = jp.Smooth(masking.antarray.vis, 0, 5, 1, Nprocess=self.Nprocess)
		#--------------------------------------------------

		# Mask manually if masktime
		# masktime => masknsplit (bool)
		gaint = self.gaint.copy()
		if (masktime is not None) : 
			masktime = jp.npfmt(masktime)
			if (masktime.size == 2) : masktime = masktime[None,:]
			# masktime is 2D, each masktime[i] is one mask range
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
			self.gaint[self.gaint.mask] = jp.ResetMasked(self.gaint, 0, self.Nprocess)
			self.gaint = self.gaint.data
		#--------------------------------------------------
		Nf = self.gaint.shape[1]
		excludelist = range(Nf/5) + range(Nf*4/5, Nf)
		self.gaintmean, self.Nlinelist_gaint, arrstd = self._Gainmean(self.gaint, 0, Nline, excludelist, gaintper, gainttimes)
		#--------------------------------------------------

		# Plot
		titlehead = r'$g(t)$'
		fighead = 'gaint'
		x = np.arange(self.gaint.shape[0])
		xlabel = 'time-axis'
		xaxes = None
		gaint = gaint[:,self.Nlinelist_gaint]
		gaint.real /= gaint.real.mean(0)[None,:,:]
		gaint.real /= self.gaintmean.real.mean(0)[None,:,:]
		try : 
			gaint.imag /= gaint.imag.mean(0)[None,:,:]
			gaint.imag /= self.gaintmean.imag.mean(0)[None,:,:]
		except : pass
		self._Plot(gaint, self.gaintmean, 0, titlehead, fighead, self.outdir, x, xlabel, xaxes, self.Nlinelist_gaint, legendsize, legendloc, progressbar)
		if (self.verbose) : print '    end   @', jp.Time(1)+'\n'





	def _Nlinelist( self, per, times, array, axis, Nline, excludelist=[] ) : 
		'''
		Nline:
			How many lines do you want to plot
			For example, plot Gaint, there are 512 freq points, total 512 lines. Set Nline=20, will just plot 20 (uniform interval) of them
			Will get Nlinelist (list of int/index)

		excludelist:
			Parameter Nline will get Nlinelist, if Nlinelist[i] in excludelist, throw it (don't plot)

		return: 
			[array, axis, Nlinelist, arrstd]
			array has been moved axis to 0
		'''
		axis, per, times = int(round(axis)), int(round(per)), int(round(times))
		if (axis not in [0,1]) : jp.Raise(Exception, 'axis must be 0=>time OR 1=>freq')
		# Move axis to 0
		array = jp.ArrayAxis(array, axis, 0, 'move')  #@#@
		#--------------------------------------------------
		# Number of curves in the figure
		try : Nline = int(round(Nline))
		except : Nline = 20
		# time or freq to be plotted
		if(Nline > array.shape[1]) : Nline = array.shape[1]
		#--------------------------------------------------
		try : excludelist = list(excludelist)
		except : excludelist = []
		#--------------------------------------------------
		# Judge which slice to be plotted
		darray = (array - jp.Smooth(array, 0, per, times)).std(0)  # real value, 2D, n1xn2
		# Reorder baseline from small to large
		ddarray = darray.std(0)
		ddarray = ddarray + 1j*np.arange(ddarray.size)
		nbl = np.sort(ddarray).imag.astype(int)[0]  # one
		darray = darray[:,nbl]
		ddarray = 0 #@
		#--------------------------------------------------
		Nlinelist = []
		nsplit=np.linspace(0, array.shape[1], Nline+1).astype(int)
		darray = darray + 1j*np.arange(darray.size) 
		for i in xrange(Nline) : 
			n = darray[nsplit[i]:nsplit[i+1]]
			n = np.sort(n).imag.astype(int)[0]
			if (n not in excludelist) : Nlinelist.append(n)
		Nlinelist = np.array(Nlinelist)
		arrstd = darray.real[Nlinelist]  # std of each line
		return [array, axis, Nlinelist, arrstd]



	def _Gainmean( self, array, axis, Nline, excludelist, per, times ) : 
		# Normalize first
		axis, per, times = int(round(axis)), int(round(per)), int(round(times))
		if (axis not in [0,1]) : jp.Raise(Exception, 'axis must be 0=>time OR 1=>freq')
		# Move axis to 0
		array = jp.ArrayAxis(array, axis, 0, 'move')  #@#@
		array.real /= array.real.mean(0)[None,:]
		try : array.imag /= array.imag.mean(0)[None,:]
		except : pass
		#--------------------------------------------------
		array, axi, Nlinelist, arrstd = self._Nlinelist(per, times, array, 0, Nline, excludelist)
		array = array[:,Nlinelist,:]  # select
		#--------------------------------------------------
		antarray = AntArray(verbose=False)
		antarray.vis = array
		masking = Masking(antarray)
		masking.MaskLoop(0, per, times, 4, 5, 0)
		array = jp.Smooth(masking.antarray.vis, 0, per, times)
		array.real /= array.real.mean(0)[None,:]
		try : array.imag /= array.imag.mean(0)[None,:]
		except : pass
		array = array.mean(1)
		#--------------------------------------------------
		if   (axis == 0) : array = array[:,None,:]
		elif (axis == 1) : array = array[None,:,:]
		return [array, Nlinelist, arrstd]





	def _Plot( self, array, arraymean, axis, titlehead, fighead, outdir, xlist=None, xlabel='', xaxes=None, Nlinelist=None, legendsize=10, legendloc=1, progressbar=None, xhour=False ) : 
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
			Used to plt_axes(xaxes[0], xaxes[1], xaxes[2])

		legendsize:
			Use for plt.legend(fontsize=legendsize)
		'''
		if (self.verbose) : 
			try : progressbar.Progress()
			except : pass
		vistype = self.antarray.vistype[:-1]
		axis = int(round(axis))
		if (axis not in [0,1]) : jp.Raise(Exception, 'axis must be 0=>time OR 1=>freq')
		# Move axis to 0
		array = jp.ArrayAxis(array, axis, 0, 'move')
		arraymean = jp.ArrayAxis(arraymean, axis, 0, 'move')
		#--------------------------------------------------
		# vmin, vmax
		vminf, vmaxf = [], []
		for i in xrange(arraymean.shape[2]) : # Plot each baseline
			ai, ami = array[:,:,i], arraymean[:,0,i]
			def _MinMax( array, ratiomax=0.5, ratiomin=0.5 ) : 
				vmin, vmax = array.min(), array.max()
				vmin, vmax = vmin-ratiomin*(vmax-vmin), vmax+ratiomax*(vmax-vmin)
				if (vmin < 0) : vmin = 0
				return vmin, vmax
			ratiomin = 0.5 if(axis==0)else 0.2
			ratiomax = 0.2
			if (vistype == 'auto') : vmin, vmax = _MinMax(ami,ratiomax, ratiomin)
			else : 
				vminr, vmaxr = _MinMax(ami.real, ratiomax, ratiomin)
				vmini, vmaxi = _MinMax(ami.imag, ratiomax, ratiomin)
				vmin, vmax = min(vminr, vmini), max(vmaxr, vmaxi)
			vminf.append(vmin)
			vmaxf.append(vmax)
		vminf, vmaxf = np.array(vminf).min(), np.array(vmaxf).max()
		#--------------------------------------------------
		for i in xrange(arraymean.shape[2]) : # Plot each baseline
			ai, ami = array[:,:,i], arraymean[:,0,i]
			if (self.verbose) : 
				try : progressbar.Progress()
				except : pass
			chanidx = self.antarray.Blorder.blorder[self.antarray.visorder][i]
			if (chanidx[0]==chanidx[1]): chanidxstr= str(chanidx[0])
			else : chanidxstr = str(chanidx[0])+'-'+str(chanidx[1])
			#--------------------------------------------------
			if (vistype == 'auto') : color = plt_color(ai.shape[1])
			else : color = plt_color(2*ai.shape[1])
			if (axis == 1) : lw = np.linspace(5, 1, ai.shape[1])
			#--------------------------------------------------
			for j in xrange(ai.shape[1]) : 
				if (vistype == 'auto') : 
					if (axis == 0) : plt.plot(xlist, ai[:,j], color=color[j], ls='', marker='o', markersize=3, label=str(Nlinelist[j]+1))
					else : plt.plot(xlist, ai[:,j], color=color[j], label=str(Nlinelist[j]+1), lw=lw[j])
				else : 
					if (axis == 0) : 
						plt.plot(xlist, ai[:,j].real, color=color[2*j], ls='', marker='o', markersize=3, label=str(Nlinelist[j]+1)+', real')
						plt.plot(xlist, ai[:,j].imag, color=color[2*j+1], ls='', marker='o', markersize=3, label=str(Nlinelist[j]+1)+', imag')
					else : 
						plt.plot(xlist, ai[:,j].real, color=color[2*j], label=str(Nlinelist[j]+1)+', real', lw=lw[j])
						plt.plot(xlist, ai[:,j].imag, color=color[2*j+1], label=str(Nlinelist[j]+1)+', imag', lw=lw[j])
			# ----- Plot arraymean -----
			lw = 3 if(axis==0)else 1
			if (vistype == 'auto') : 
				if   (ai.shape[1] == 1) : color = 'r'
				elif (ai.shape[1] == 2) : color = 'g'
				elif (ai.shape[1] == 3) : color = 'k'
				elif (ai.shape[1] == 4) : color = 'c'
				elif (ai.shape[1] == 5) : color = 'y'
				else : color = 'b'
				plt.plot(xlist, ami, color=color, lw=lw, label='mean')
			else : 
				if   (ai.shape[1] == 1) : color = ['c', 'g']
				elif (ai.shape[1] == 2) : color = ['c', 'y']
				else : color = ['c', 'b']
				plt.plot(xlist, ami.real, color=color[0], lw=lw, label='mean, real')
				plt.plot(xlist, ami.imag, color=color[1], lw=lw, label='mean, imag')
			plt.legend(fontsize=legendsize, loc=legendloc)
			plt.xlabel(xlabel, size=16)
			plt.xlim(int(round(xlist.min())),int(round(xlist.max())))
			try : plt_axes(xaxes[0], xaxes[1], xaxes[2])
			except : pass
			plt.ylabel('[A.U.]', size=16)
			plt.ylim(vminf, vmaxf)
			if (vistype == 'auto') : 
				nstr, channame = str(ai.shape[1])+'+1', 'channel='
			else : 
				nstr, channame = str(2*ai.shape[1])+'+2', 'baseline='
			plt.title(titlehead+' of '+channame+chanidxstr+', total '+nstr+' curves', size=16)
			plt.savefig(outdir+fighead+'_chan'+chanidxstr+'_N'+str(ai.shape[1])+'.png')
			plt.close()





	def Save( self, outname=None ) : 
		''' outname: Absolute path of output .hdf5 '''
		if (outname is not None) : self.outname = str(outname)
		if (self.outname[-5:] !='.hdf5') : self.outname += '.hdf5'
		if (self.verbose) : 
			print 'CaliGain.Save'
			print '    Saving to  "'+jp.AbsPath(self.outname)+'"\n'
		caligain = CaliGain(verbose=False)
		caligain.__dict__ = self.__dict__
		caligain.__dict__.pop('antarray')
		caligain.__dict__.pop('masking')
		classhdf5 = jp.ClassHdf5(caligain, self.outname, verbose=False)
		classhdf5.Save()





	def Read( self, outname=None ) : 
		''' outname: Absolute path of output .hdf5 '''
		if (outname is not None) : self.outname = str(outname)
		if (self.outname[-5:] !='.hdf5') : self.outname += '.hdf5'
		if (self.verbose) : 
			print 'CaliGain.Read'
			print '    Reading from  "'+jp.AbsPath(self.outname)+'"\n'
		caligain = CaliGain(verbose=False)
		classhdf5 = jp.ClassHdf5(caligain, self.outname, verbose=False)
		classhdf5.Read()
		caligain.__dict__.pop('antarray') 
		caligain.__dict__.pop('masking') 
		self.__dict__.update(caligain.__dict__)


