import jizhipy as jp
from jizhipy.Plot import *
from AntArray import *
from Masking import *
##################################################



def PhaseBrightsource( x, y, lat, Dec, freq, RA0, RA ) : 
	'''
	NOTE THAT:
		x, y, lat, Dec, freq, RA0 must be ONE
		RA can be ONE or N-D Array

	     North             +y
	West       East     -x    +x     zenoth: +z
        South             -y

	x, y:
		Coordinate (above) of the baseline

	lat: 
		in rad
		Latitude of the antenna array/site

	Dec:
		in rad
		Declination of the antenna pointing (transit)

	freq:
		in MHz

	RA0:
		in rad
		Reference of the RA 
		For bright source calibration, RA0 = RAsource

	RA:
		in rad
		Transit RA
		Can be one or N-D array
	'''
	Xb = x
	Yb = y *np.sin(lat)
	Zb = y *np.cos(lat)
	#--------------------
	Xs =  np.cos(Dec) *np.sin(RA-RA0)
	Ys = -np.cos(Dec) *np.cos(RA-RA0)
	Zs =  np.sin(Dec)
	#--------------------
	phase = 2*np.pi/(300./freq) * (Xb*Xs + Yb*Ys + Zb*Zs)
	return phase





def _DoMultiprocess_FitBeam( iterable ) : 
	n1, n2  = iterable[0][0]
	beam = iterable[1]
	timem, s0, Dec, verbose = iterable[2]
	if (verbose) : progressbar = jp.ProgressBar('    pid='+str(os.getpid())+':', len(beam), False)
	#--------------------------------------------------
	def func(x, p) :  # x is time, not angle
		beta=jp.Circ2Sph((x-p[1])/60.*15*np.pi/180, Dec*np.pi/180)
		y = p[0] * np.exp(-beta**2 /2 /p[2]**2) +p[3]
	#	if (p[3] < 0) : y = -1e10
		return y
	#--------------------------------------------------
	pf = []
	for i in xrange(n2-n1) : 
		if (verbose) : progressbar.Progress()
		p0 = [beam[i].max(), timem[beam[i]==beam[i].max()].mean(), s0, 0]
		p = jp.FuncFit(func, timem, beam[i], p0, warning=False)[0]
		pf.append( p )
	return jp.npfmt(pf)



def _DoMultiprocess_FitPhase( iterable ) : 
	n1, n2 = iterable[0][0]
	phase, t0list = iterable[1]
	timem, Dec, freq, bl, verbose, Nf, lat, fitLew = iterable[2]
	if (verbose) : progressbar = jp.ProgressBar('    pid='+str(os.getpid())+':', len(phase), False)
	lat, Dec = lat*np.pi/180, Dec*np.pi/180
	RA  =  timem/60.*15*np.pi/180
	RA0 = t0list/60.*15*np.pi/180
	#--------------------------------------------------

	pf = []
	for i in xrange(n2-n1) : 
		if (verbose) : progressbar.Progress()
		nv, nf = (i+n1)/Nf, (i+n1)%Nf
		x, y = bl[nv,:2]  #@#@
		#--------------------------------------------------
		def funcNO( phi, p ) : 
			phase = PhaseBrightsource(x, y, lat, Dec, freq, RA0, phi) + p[0]
			return np.sin(phase)
		#--------------------------------------------------
		def funcYES( phi, p ) : 
			phase = PhaseBrightsource(p[0], y, lat, Dec, freq, RA0, phi) + pf[-1][0]
			return np.sin(phase)
		#--------------------------------------------------
		def funcBOTH( phi, p ) : 
			phase = PhaseBrightsource(p[1], y, lat, Dec, freq, RA0, phi) + p[0]
			return np.sin(phase)
		#--------------------------------------------------

		pf.append( list(jp.FuncFit(funcNO, RA, np.sin(phase[i]), [1])[0] %(2*np.pi)) )
		#--------------------------------------------------

		if (fitLew == 0) : pf[-1] += [x]
		elif (fitLew == 1) : 
			pf[-1] += list(jp.FuncFit(funcYES, RA, np.sin(phase[i]), [x])[0])
		elif (fitLew == 2) : 
			pf[-1] = jp.FuncFit(funcBOTH, RA, np.sin(phase[i]), [pf[-1][0], x])[0]
		#--------------------------------------------------
	return jp.npfmt(pf)





def _DoMultiprocess_FitPhase( iterable ) : 
	n1, n2 = iterable[0][0]
	phase, t0list = iterable[1]
	timem, Dec, freq, bl, verbose, Nf, lat, fitLew = iterable[2]
	if (verbose) : progressbar = jp.ProgressBar('    pid='+str(os.getpid())+':', len(phase), False)
	#--------------------------------------------------
	pf = []
	for i in xrange(n2-n1) : 
		t0 = t0list[i] #@
		if (verbose) : progressbar.Progress()
		nv, nf = (i+n1)/Nf, (i+n1)%Nf
		#--------------------------------------------------
		def funcNO(x, p) : 
			beta=jp.Circ2Sph((x-t0)/60.*15*np.pi/180,Dec*np.pi/180)
			angle= 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
		#	angle= -2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0] # (2)
		#	angle= 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) + bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0] # (4)
		#	angle= -2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) + bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0] # (1)
			return np.sin(angle)
		#--------------------------------------------------
		def funcYES(x, p) : 
			beta=jp.Circ2Sph((x-t0)/60.*15*np.pi/180,Dec*np.pi/180)
			angle = 2*np.pi/(300./freq[nf]) *( p[0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + pf[-1][0]
			return np.sin(angle)
		#--------------------------------------------------
		def funcBOTH(x, p) : 
			beta = jp.Circ2Sph((x-t0)/60.*15*np.pi/180, Dec*np.pi/180)
			angle = 2*np.pi/(300./freq[nf]) *( p[1]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
			return np.sin(angle)
		#--------------------------------------------------
		pf.append( list(jp.FuncFit(funcNO, timem, np.sin(phase[i]), [1])[0] %(2*np.pi) ) )
		#--------------------------------------------------
		if (fitLew == 0) : pf[-1] += [bl[nv,0]]
		elif (fitLew == 1) : 
			pf[-1] += list(jp.FuncFit(funcYES, timem, np.sin(phase[i]), [bl[nv,0]])[0])
		elif (fitLew == 2) : 
			pf[-1] = jp.FuncFit(funcBOTH, timem, np.sin(phase[i]), [pf[-1][0], bl[nv,0]])[0]
		#--------------------------------------------------
	return jp.npfmt(pf)





def _DoMultiprocess_FitVis( iterable ) : 
	n1, n2 = iterable[0][0]
	vis, A0, t0list = iterable[1]
	x0, s0, Dec, freq, bl, showprogress, Nv, lat = iterable[2]
	if (showprogress) : progressbar = jp.ProgressBar('CaliPhase.FitPhase(), pid='+str(os.getpid())+':', len(phase))
	pf = []
	#--------------------------------------------------
	for i in xrange(n2-n1) : 
		if (showprogress) : progressbar.Progress()
		nf, nv = (i+n1)/Nv, (i+n1)%Nv
		def func(x, p) : 
			beta = jp.Circ2Sph((x-t0list[i])/3600.*15*np.pi/180, Dec*np.pi/180)
			A = p[0] * np.exp(-beta**2 /2 /p[1]**2)
			angle = 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[2]
			return A*np.sin(angle)
		#	return A*np.cos(angle)
		#--------------------------------------------------
		p0 = [A0[i], s0, 1]
		pf.append( jp.FuncFit(func, x0, vis[i], p0)[0] )
	return jp.npfmt(pf)  # shape=(n2-n1, 4)


def _DoMultiprocess_FitExp( iterable ) : 
	n1, n2 = iterable[0][0]
	vis, t0list = iterable[1]
	x0, Dec, freq, bl, showprogress, Nv, lat = iterable[2]
	if (showprogress) : progressbar = jp.ProgressBar('CaliPhase.FitExp(), pid='+str(os.getpid())+':', len(phase))
	pf = []
	#--------------------------------------------------
	for i in xrange(n2-n1) : 
		if (showprogress) : progressbar.Progress()
		nf, nv = (i+n1)/Nv, (i+n1)%Nv
		beta = jp.Circ2Sph((x0-t0list[i])/3600.*15*np.pi/180, Dec*np.pi/180)
		A = abs(vis[i])
		def func(x, p) : 
			angle = 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
			return A*np.sin(angle)
		#	return A*np.cos(angle)
		#--------------------------------------------------
		p0 = np.pi
		pf.append( jp.FuncFit(func, x0, vis.imag[i], p0)[0] )
	return jp.npfmt(pf).flatten()



def _DoMultiprocess_Plot( iterable ) : 
	n1, n2 = iterable[0][0]
	array = iterable[1][:-1]
	strbl = iterable[1][-1]
	freq, func, color, outdir, which, nsigma, freqminmax = iterable[2][:7]
	freqminmax = jp.npfmt(freqminmax).round().astype(int)
	which = which.lower()
	if (which == 'deff') : dyDeff, fwhm2sigmafactor, antform = iterable[2][-3:]
	elif (which == 'lew') : fitLew, bl, dyLew = iterable[2][-3:]
	elif (which in ['phaseadd', 'phase']) : fitLew, bl = iterable[2][-2:]
	#--------------------------------------------------
	for i in xrange(n2-n1) : # baseline
		if (which == 'amp') : arraymax = np.array([0.,0])+1e30
		elif (which in ['lew', 'phaseadd', 'phase']) : strl = '=(%.3f, %.3f)' % tuple(bl[i][:2])
		amean = []
		for j in xrange(len(array)) : # Different method: FitBeam, FitVis, ...
			if (array[j] is None) : continue
			nmf = array[j][i].size/5
			nmf = nmf if(nmf%2==1) else nmf+1
			if (which == 'amp') : 
				arrayfilt = spsn.medfilt(array[j][i], nmf)
				arraymax[j] = arrayfilt.max()
				amean.append(arrayfilt.mean())
			else : 
				n1, n2 = len(freq)/5, len(freq)*4/5
				amean.append(array[j][i][n1:n2].mean())
			plt.plot(freq, array[j][i], color=color[j], ls='', marker='o', markersize=3, label=func[j])
		amean = np.array(amean).mean()
		plt.xlim(freqminmax[0], freqminmax[1])
		plt.xlabel(r'$\nu$ [MHz]', size=16)
		plt_axes('x', 'both', [25,5])
		#--------------------------------------------------
		if (which == 'amp') : 
			plt.ylabel('Amp [A.U.]', size=16)
			plt.ylim(0, 1.4*arraymax.min())
		#	plt_axisformat('y')
			plt.legend()
			plt.title('Fitted amplitude of abs(visibility), baseline='+strbl[i], size=16)
			plt.savefig(outdir+'Amp_'+strbl[i]+'_'+str(nsigma)+'.png')
		#--------------------------------------------------
		elif (which == 'deff') : 
			plt.ylabel(r'Deff$=\frac{\lambda}{%.3f \cdot \sigma}$ [m]' % fwhm2sigmafactor, size=16)
			if (dyDeff is not None) : plt.ylim(amean-dyDeff, amean+dyDeff)
			plt_axes('y', 'both', [0.1, 0.05], '%.1f')
			plt.legend()
		#	plt.title(r'Fitted effective diameter Deff$=\frac{\lambda}{'+('%.3f' % fwhm2sigmafactor)+r' \cdot \sigma}$, baseline='+strbl[i], size=16)
			plt.title(r'Fitted effective diameter ('+str(antform)+'m), baseline='+strbl[i], size=16)
			plt.savefig(outdir+'Deff_'+strbl[i]+'_'+str(nsigma)+'.png')
		#--------------------------------------------------
		elif (which == 'lew') : 
			plt.ylabel(r'$L_{ew}$ [m]', size=16)
			if (dyLew is not None) : plt.ylim(amean-dyLew, amean+dyLew)
			plt_axes('y', 'both', [0.05, 0.01], '%.2f')
			plt.legend()
			plt.title(r'Fitted $L_{ew}$, baseline='+strbl[i]+strl,size=16)
			plt.savefig(outdir+'Lew_'+strbl[i]+'_'+str(nsigma)+'_m'+str(fitLew)+'.png')
		#--------------------------------------------------
		elif (which == 'phaseadd') : 
			plt.ylabel('Phaseadd [deg]', size=16)
			plt_axes('y', 'both', [30,5])
			plt.ylim(0, 360)
			plt.title('Fitted additional phase, baseline='+strbl[i]+strl, size=16)
			plt.legend()
			plt.savefig(outdir+'Phaseadd_'+strbl[i]+'_'+str(nsigma)+'_m'+str(fitLew)+'.png')
		plt.close()


##################################################
##################################################
##################################################





class CaliPhase( object ) : 



	def __init__( self, caligain, fwhm2sigmafactor=2.287, outdir='' ) : 
		'''
		antarray: instance of class:AntArray ! np.ndarray
		masking:  instance of class:Masking
		caligain: instance of class:CaliGain

		fwhm2sigmafactor:
			sigma = lambda /Deff /fwhm2sigmafactor
			Gaussian beam: 2.287
			Reza: 2
		'''
		self.Nprocess, self.verbose = caligain.Nprocess, caligain.verbose
		if (self.verbose) : print '-------------------- CaliPhase --------------------\n'
		if (fwhm2sigmafactor is not None) : 
			try : self.fwhm2sigmafactor = float(fwhm2sigmafactor)
			except : self.fwhm2sigmafactor = 2.287
		#--------------------------------------------------
		self.caligain = caligain
		self.masking  = self.caligain.masking
		self.antarray = self.masking.antarray
		self.vistype  = self.masking.vistype[:-1]
		try : a, b = self.antarray.File.transittime, self.antarray.File.nfiletransit
		except: jp.Raise(Exception,"NOT exist fo['transitsource'], you can set it by hand in AntArray.__init__(transittime=)\n")
		#--------------------------------------------------
		self.outdir = jp.Outdir((None,'file'), (0,'file'))
		try : self.outdir = jp.Mkdir(self.outdir + outdir)
		except : pass
		self.timeper, self.timetimes = caligain.timeper, caligain.timetimes





	def WhichTransitsource( self, which ) : 
		''' which: int, 0,1,2,... '''
		which = jp.Num(which)
		try : 
			self.RA, self.Dec = self.antarray.File.transitRADec[which].flatten()
			self.whichtransitsource = which
		except : jp.Raise(Exception, "self.antarray.File doesn't have .transitRADec OR self.antarray is np.ndarray. Please use CaliPhase.PointDec()")



#	def PointDec()





	def Fringe( self, nsigma=4, plotfreq=None, plotnv=None, nyshift=2, yshiftpos=None ) : 
		'''
		nyshift:
			self.yshift = self.vis(10sigma)[:nyshift*1sigma]
			nyshift = 0: don't shift

		yshiftpos:
			False, None, 'left', 'right', 'both'
			False: don't shift
			None: use existed self.yshift
			'left' : self.yshift=self.vis(10sigma)[:nyshift*1sigma]
			'right': self.yshift=self.vis(10sigma)[-nyshift*1sigma:]
			'both' : (self.yshift=self.vis(10sigma)[:nyshift*1sigma]+self.yshift=self.vis(10sigma)[-nyshift*1sigma:])/2

		Take where is the fringe of the bright source
		Generate: 
			self.timerange: pixels range of the fringe
			self.vis: MaskedArray of fringe data, has been reset by maskvalue
			self.timem: in minute, self.timem.size=self.vis.shape[0]
			self.yshift: a constant shift of y-axis


		If there is an extra constant input/interference in front of A/D, this constant will become delta function (real value) after FFT, that will shift the real part of visibilities by a constant, but doesn't affect the imaginary part.

		For auto-correlation, we can't remove this shifting because we don't know how much it is, however, for cross-correlation, the base line/average? should be zero because of exp{ix}, we can remove the shifting in the cross-correlation.
		'''
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliPhase.Fringe: start @', starttime
		try : self.whichtransitsource+0
		except : 
			if (len(self.antarray.File.transitname) == 1) : self.WhichTransitsource(0)
			else : jp.Raise(Exception, 'there are '+str(len(self.antarray.File.transitname))+' transitsources='+str(self.antarray.File.transitname)+', but now self.whichtransitsource=None, please select one of them using CaliPhase.WhichTransitsource()')
		if (nsigma is not None) : 
			try : nsigma = float(nsigma)
			except : nsigma = 4
			if (nsigma == int(nsigma)) : nsigma = int(nsigma)
		else : nsigma = 4
		self.nsigma = 10
		self.plotfreq, self.plotnv = plotfreq, plotnv
		inttime = self.antarray.Ant.inttime
		freq = self.antarray.Ant.freq[self.antarray.freqorder]
		antform = jp.npfmt(self.antarray.Ant.antform)[0]
		self.vis = np.ma.MaskedArray(self.antarray.vis, self.masking.mask)
		#--------------------------------------------------

		pixsource = int(self.antarray.File.transittime[self.whichtransitsource] / inttime) - self.antarray.File.nfile[0] * self.antarray.Ant.N0
		#--------------------------------------------------

		sigma = 300/freq.mean()/antform/self.fwhm2sigmafactor
		sigma = jp.Sph2Circ(sigma, self.Dec*np.pi/180)  # rad
		# Total pixels of the fringe in 1 sigma
		pix1 = int(round(24*3600/(2*np.pi)*sigma/inttime)) # int

		pix10 = [pixsource-10*pix1/2, pixsource+10*pix1/2] #(n1,n2)
		if (pix10[0] < 0) : pix10[0] = 0
		if (pix10[1] > self.antarray.vis.shape[0]) : pix10[1] = self.antarray.vis.shape[0]
		self.pix10sigma = pix10
		self.timem = self.antarray.timem[pix10[0]:pix10[1]]
		if (self.vistype == 'cross') : 
			if (nyshift is 0 or nyshift is None) : yshiftpos = None
			if (yshiftpos is False) : self.yshift = np.zeros((1,)+self.vis.shape[1:], complex)
			elif (yshiftpos is None) : 
				try : self.yshift + 0
				except : self.yshift = np.zeros((1,)+self.vis.shape[1:], complex)
			elif (yshiftpos == 'left') : 
				self.yshift = self.vis[pix10[0]:pix10[1]][:nyshift*pix1].mean(0)[None,:,:]
			elif (yshiftpos == 'right') : 
				self.yshift = self.vis[pix10[0]:pix10[1]][-nyshift*pix1:].mean(0)[None,:,:]
			else : 
				self.yshift = np.append(self.vis[pix10[0]:pix10[1]][:nyshift*pix1], self.vis[pix10[0]:pix10[1]][-nyshift*pix1:], 0).mean(0)[None,:,:]
			self.vis = self.vis - self.yshift
			self.antarray.vis -= self.yshift
		#--------------------------------------------------

		if (self.vistype == 'cross') : self.PlotFringe(None, None, self.vis[pix10[0]:pix10[1]].data, ['Real', 'Imaginary', 'Abs'], ['b', 'r', 'g'])
		else : self.PlotFringe(None, None, self.vis[pix10[0]:pix10[1]].data, ['Auto'], ['b'])
		#--------------------------------------------------
		
		self.nsigma = nsigma
		pixf = [pixsource-self.nsigma*pix1/2, pixsource+self.nsigma*pix1/2] # (n1, n2)
		if (pixf[0] < 0) : pixf[0] = 0
		if (pixf[1] > self.antarray.vis.shape[0]) : pixf[1] = self.antarray.vis.shape[0]
		self.vis = self.vis[pixf[0]:pixf[1]]
		self.timem = self.antarray.timem[pixf[0]:pixf[1]]
		#--------------------------------------------------

		if (self.vistype == 'cross') : self.PlotFringe(None, None, self.vis.data, ['Real', 'Imaginary', 'Abs'], ['b', 'r', 'g'])
		else : self.PlotFringe(None, None, self.vis.data, ['Auto'], ['b'])
		#--------------------------------------------------

		if (self.verbose) : 
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliPhase.Fringe:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'




	def Smooth( self, timeper=None, timetimes=None, freqper=None, freqtimes=None ) : 
		'''
		jp.Smooth(self.vis, 1, freqper, freqtimes)
		jp.Smooth(self.vis, 0, timeper, timetimes)
		'''
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliPhase.Smooth: start @', starttime
		if (timeper is None) : timeper = self.timeper
		if (timetimes is None) : timetimes = self.timetimes
		self.timeper, self.timetimes = timeper, timetimes
		if (freqper is None) : freqper = 4
		if (freqtimes is None) : freqtimes = 1
		self.freqper, self.freqtimes = freqper, freqtimes
		#--------------------------------------------------
		if (timeper % 2 == 1) : n = (timeper-1) * timetimes /20
		else : n = timeper * timetimes /20
		array = [self.vis[n:-n]]
		#--------------------------------------------------
		mask = self.vis.mask.copy()
		self.vis = jp.Smooth(self.vis.data, 1, freqper, freqtimes, Nprocess=self.Nprocess)
		self.vis = jp.Smooth(self.vis, 0, timeper, timetimes, Nprocess=self.Nprocess)
		self.vis = np.ma.MaskedArray(self.vis, mask)
		self.vis = self.vis[n:-n]
		self.timem = self.timem[n:-n]
		#--------------------------------------------------
		array.append(self.vis)
		if (self.vistype == 'cross') : self.PlotFringe(None, None, array, [['Real','Imaginary','Abs'],['Smoothed','Smoothed','Smoothed']], [['b','r','k'],['m','g','c']], [1,3])
		else : self.PlotFringe(None, None, array, [['Auto-corr'],['Smoothed']], [['b'],['r']], [1,3])
		if (self.verbose) : 
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliPhase.Smooth:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'



	def FitBeam( self ) : 
		''' Fit the beam with abs(vis) '''
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliPhase.FitBeam: start @', starttime
			if (self.vistype == 'auto') : 
				print '    Warning: CaliPhase.vistype="'+self.vistype+'", NOT "cross"'
				print '    Do nothing !'
				endtime = jp.Time(1)
				costtime = jp.Time(starttime, endtime)
				tottime = jp.Time(self.starttime, endtime)
				print 'CaliPhase.FitBeam:  end  @', endtime
				print 'Cost:', costtime+'     Total:', tottime+'\n'
				return
			else : print '    NOTE THAT: when using CaliPhase.FitBeam, nsigma should be enough large to cover the whole fringe/beam profile'
		# Initial guess
		inttime = self.antarray.Ant.inttime
		antform = jp.npfmt(self.antarray.Ant.antform)[0]
		s0 = 300/self.antarray.Ant.freq.mean()/antform/2.287
		#--------------------------------------------------
		beam = abs(self.vis.data)
	#	beam = ArrayAxis(beam, 0, -1, 'move')
		beam = beam.T  # Faster
		shape = beam.shape
		# Reshape to 2D
		beam = beam.reshape(np.prod(shape[:2]), shape[-1])
		#--------------------------------------------------
		bcast = (self.timem, s0, self.Dec, self.verbose)
		pool = jp.PoolFor(0, len(beam), self.Nprocess)
		pf = pool.map_async(_DoMultiprocess_FitBeam, beam, bcast)
		print
		#--------------------------------------------------
		pf = np.concatenate(pf)  # (nv*nf, 4)
		pf = pf.reshape(shape[:2]+(pf.shape[-1],)).T  # (4,nf,nv)
		self.Amp1b, self.Timeb, self.Sigmab, self.Amp2b = pf
		self.Ampb = self.Amp1b + self.Amp2b
		# FitBeam: Beam = Amp1b * exp(-x**2/2/s**2) + Amp2b
		# But in fact, for FitPhase, Amp=Amp1+Amp2, vis=(Amp1+Amp2)*exp(i*...)
		del pool, beam
		self.Sigmab = abs(self.Sigmab)
		#--------------------------------------------------
		np.save(self.outdir+'Amp1b_'+str(self.nsigma)+'.npy', self.Amp1b)
		np.save(self.outdir+'Amp2b_'+str(self.nsigma)+'.npy', self.Amp2b)
		np.save(self.outdir+'Ampb_'+str(self.nsigma)+'.npy', self.Ampb)
		np.save(self.outdir+'Timeb_'+str(self.nsigma)+'.npy',self.Timeb)
		np.save(self.outdir+'Sigmab_'+str(self.nsigma)+'.npy', self.Sigmab)
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		self.Phasens = 2*np.pi/300*self.antarray.Ant.freq[self.antarray.freqorder][:,None]*(-bl[:,1][None,:])*np.sin((self.Dec-self.antarray.Ant.lonlat[1])*np.pi/180)
		np.save(self.outdir+'Phasens.npy', self.Phasens)
		# vis * exp(-1j*Phasens) * exp(-1j*Phaseadd)
		#--------------------------------------------------
		# Plot
		nf, strfreq, nv, strblo, strbli, blname = self.antarray.Plotnfnv(None,None,False)
		def func(x, p) :  # x is time, not angle
			beta = jp.Circ2Sph((x-p[1])/60.*15*np.pi/180, self.Dec*np.pi/180)
			y = p[0] * np.exp(-beta**2 /2 /p[2]**2) +p[3]
			return y
		if (self.timem[-1]-self.timem[0] <= 60) : xmajor = 5
		else : xmajor = 10
		p = [self.Amp1b[nf,nv], self.Timeb[nf,nv], self.Sigmab[nf,nv], self.Amp2b[nf,nv]]
		beam = abs(self.vis.data[:,nf,nv])
		plt.plot(self.timem, beam, 'b-', label='Data')
		plt.plot(self.timem, func(self.timem, p), 'r-', label='Fitting')
		plt.legend()
		plt.xlabel('time [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [xmajor,1])
	#	plt.ylabel('[A.U.]', size=16)
		plt.title('CaliPhase.FitBeam, '+strfreq+', '+blname+strblo+strbli, size=14)
		plt.savefig(self.outdir+'CaliPhase.FitBeam_'+strfreq+'_'+strblo+'.png')
		plt.close()
		#--------------------------------------------------
		if (self.verbose) : 
			print '    Amp1b_'+str(self.nsigma)+'.npy, Amp2b_'+str(self.nsigma)+'.npy, Ampb_'+str(self.nsigma)+'.npy, Timeb_'+str(self.nsigma)+'.npy, Sigmab_'+str(self.nsigma)+'.npy, Phasens_.npy  saved to '+self.outdir
			print '    Plotting CaliPhase.FitBeam_'+strfreq+'_'+strblo+'.png'
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliPhase.FitBeam:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'





	def FitPhase( self, fitLew=0 ) : 
		'''
		(1) Because there is dPhase, fringe is not sensitive enough to Lew, we need to first fit the initial guess of dPhase.
		(2) Also because fringe is not sensitive enough to Lew, different range of fringe (nsigma=1, 2, 3, ...) have very different Lew.

		fitLew:
			=0, 1, 2. Different metohds to fit Lew
			==0: don't fit Lew, just fit Phaseadd
			==1: fit Lew and Phaseadd separately
			==2: fit Lew and Phaseadd together (not good)
		'''
		try : self.fitLew = int(fitLew)
		except : self.fitLew = 0
		#--------------------------------------------------
		if ('Timeb' not in self.__dict__.keys()) : self.FitBeam()
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliPhase.FitPhase: start @', starttime
			if (self.vistype == 'auto') : 
				print '    Warning: CaliPhase.vistype="'+self.vistype+'", NOT "cross"'
				print '    Do nothing !'
				endtime = jp.Time(1)
				costtime = jp.Time(starttime, endtime)
				tottime = jp.Time(self.starttime, endtime)
				print 'CaliPhase.FitPhase:  end  @', endtime
				print 'Cost:', costtime+'     Total:', tottime+'\n'
				return
			else : print '    NOTE THAT: when using CaliPhase.FitPhase, nsigma should be smaller to avoid the interference from other sources'
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		freq = self.antarray.Ant.freq[self.antarray.freqorder]
		inttime = self.antarray.Ant.inttime
		#--------------------------------------------------
		phase = np.angle(self.vis.data)
	#	phase = ArrayAxis(phase, 0, -1, 'move')
		phase = phase.T
		shape = phase.shape
		# Reshape to 2D
		phase = phase.reshape(np.prod(shape[:2]), shape[-1])
		#--------------------------------------------------
		pool = jp.PoolFor(0, len(phase), self.Nprocess)
		send = (phase, self.Timeb.T.flatten())
		bcast = (self.timem, self.Dec, freq, bl, self.verbose, shape[1], self.antarray.Ant.lonlat[1], self.fitLew)
		pf = pool.map_async(_DoMultiprocess_FitPhase, send, bcast)
		print
		#--------------------------------------------------
		pf = np.concatenate(pf)  # (nv*nf, 2)
		pf = pf.reshape(shape[:2]+(pf.shape[-1],)).T  # (2,nf,nv)
		self.__dict__['Phaseaddp'+str(self.fitLew)] = pf[0] %(2*np.pi)
		self.__dict__['Lewp'+str(self.fitLew)] = pf[1]
		#--------------------------------------------------
		np.save(self.outdir+'Lewp'+str(self.fitLew)+'_'+str(self.nsigma)+'.npy', self.__dict__['Lewp'+str(self.fitLew)])
		np.save(self.outdir+'Phaseaddp'+str(self.fitLew)+'_'+str(self.nsigma)+'.npy', self.__dict__['Phaseaddp'+str(self.fitLew)])
		#--------------------------------------------------
		# Plot
		nf, strfreq, nv, strblo, strbli, blname = self.antarray.Plotnfnv(None,None,False)
		def func(x, p) : 
			beta = jp.Circ2Sph((x-p[1])/60.*15*np.pi/180, self.Dec*np.pi/180)
			A = p[0] * np.exp(-beta**2 /2 /p[2]**2)
			angle = 2*np.pi/(300./freq[nf]) *( p[4]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((self.Dec-self.antarray.Ant.lonlat[1])*np.pi/180) ) + p[3]
			return A*np.sin(angle)
		if (self.timem[-1]-self.timem[0] <= 60) : xmajor = 5
		else : xmajor = 10
		p = (self.Ampb[nf,nv], self.Timeb[nf,nv], self.Sigmab[nf,nv], self.__dict__['Phaseaddp'+str(self.fitLew)][nf,nv], self.__dict__['Lewp'+str(self.fitLew)][nf,nv])
		vis = self.vis.data[:,nf,nv]
		plt.plot(self.timem, vis.imag, 'b-', label='vis.imag')
		plt.plot(self.timem, func(self.timem, p), 'r-', label='Fitting, m'+str(self.fitLew))
		plt.legend()
		plt.xlabel('time [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [xmajor,1])
	#	plt.ylabel('[A.U.]', size=16)
		plt.ylim(-self.Ampb[nf,nv]*1.1, self.Ampb[nf,nv]*1.1)
		plt.title('CaliPhase.FitPhase, '+strfreq+', '+blname+strblo+strbli, size=14)
		plt.savefig(self.outdir+'CaliPhase.FitPhase_m'+str(fitLew)+'_'+strfreq+'_'+strblo+'.png')
		plt.close()
		#--------------------------------------------------
		if (self.verbose) : 
			print '    Phaseaddp'+str(self.fitLew)+'_'+str(self.nsigma)+'.npy, Lewp'+str(self.fitLew)+'_'+str(self.nsigma)+'.npy  saved to '+self.outdir
			print '    Plotting CaliPhase.FitPhase_m'+str(fitLew)+'_'+strfreq+'_'+strblo+'.png'
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliPhase.FitPhase:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'





#	def FitVis(self, Nprocess=None): 
#		if ('Timeb' not in self.__dict__.keys()) : self.FitBeam()
#		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
#		freq = self.antarray.Ant.freq
#		inttime = self.antarray.Ant.inttime
#		s0 = 300/freq.mean()/self.antarray.Ant.dishdiam/2.287
#		A0 = abs(self.vis.data).max(0).flatten()
#		x0 = np.arange(self.vis.shape[0])*inttime  # second
#		#--------------------------------------------------
#		# Move time axis to -1
#	#	vis = ArrayAxis(self.vis.data.imag, 0, -1, 'move')
#		vis = vis.data.imag.T
#		shape = vis.shape  # (512,528,3600)
#		# Reshape to 2D
#		vis = vis.reshape(np.prod(shape[:2]), shape[-1])
#		#--------------------------------------------------
#		pool = jp.PoolFor(0, len(vis), Nprocess)
#		send = (vis, A0, self.Timeb.flatten())
#		bcast = (x0, s0, self.Dec, freq, bl, self.showprogress, shape[1], self.antarray.Ant.lonlat[1])
#		pf = pool.map_async(_DoMultiprocess_FitVis, send, bcast)
#		#--------------------------------------------------
#		pf = np.concatenate(pf).T  # shape=(4, Nf*Nv)
#		del pool, bcast, vis, A0
#		pf = pf.reshape( (len(pf),) + shape[:2] )
#		self.Ampv, self.Sigmav, self.Phaseaddv = pf
#		self.Sigmav = abs(self.Sigmav)
#		self.Phaseaddv[self.Ampv<0] += np.pi
#		self.Phaseaddv %= (2*np.pi)
#		self.Ampv = abs(self.Ampv)
#		#--------------------------------------------------
#		np.save(self.outdir+'Ampv_'+str(self.nsigma)+'.npy', self.Ampv)
#		np.save(self.outdir+'Sigmav_'+str(self.nsigma)+'.npy', self.Sigmav)
#		np.save(self.outdir+'Phaseaddv_'+str(self.nsigma)+'.npy', self.Phaseaddv)
#		#--------------------------------------------------
#		# Plot
#		if (self.plotnv is None) : 
#			bl = abs(bl[:,0]) # longest East-West baseline
#			nv = np.where(bl==bl.max())[0][0]
#		else : nv = self.plotnv
#		nvbl = np.arange(len(self.antarray.Blorder.blorder))[self.antarray.visorder][nv]
#		bl = self.antarray.Blorder.Order2Bl(nvbl)
#		strbl = str(bl[0])+'-'+str(bl[1])
#		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
#		if (plotfreq is None) : plotfreq = freq.mean()
#		nf = abs(freq-plotfreq)
#		freqstr = str(int(round(plotfreq)))
#		nf = np.where(nf==nf.min())[0][0]
#		vis = self.vis[:,nf,nv]  # real+imag+masked
#		vmax = abs(vis).max()*1.05
#		#--------------------
#		def func(x, p) : 
#			beta = jp.Circ2Sph((x-p[1])/3600.*15*np.pi/180, self.Dec*np.pi/180)
#			A = p[0] * np.exp(-beta**2 /2 /p[2]**2)
#			angle = 2*np.pi/(300./freq[nf]) * ( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((self.Dec-self.antarray.Ant.lonlat[1])*np.pi/180) ) + p[3]
#			return A*np.exp(1j*angle)  # complex
#		#--------------------
#		visf = func(x0, [self.Ampv[nf,nv], self.Timeb[nf,nv], self.Sigmav[nf,nv], self.Phaseaddv[nf,nv]])
#		plt.figure(figsize=(17,6))
#		plt.subplot(1,2,1)
#		plt.plot(self.timem, vis.real, 'b-', label='data')
#		plt.plot(self.timem, visf.real, 'r-', label='Fitting')
#		plt.legend(fontsize=16)
#		plt.xlabel('$t$ [min]', size=16)
#		plt.xlim(self.timem.min(), self.timem.max())
#		plt_axes('x', 'both', [10,1])
#		plt.ylabel('[A.U.]', size=16)
#		plt.ylim(-vmax, vmax)
#		plt.title('Real part', size=16)
#		plt.subplot(1,2,2)
#		plt.plot(self.timem, vis.imag, 'b-', label='data')
#		plt.plot(self.timem, visf.imag, 'r-', label='Fitting')
#		plt.legend(fontsize=16)
#		plt.xlabel('$t$ [min]', size=16)
#		plt.xlim(self.timem.min(), self.timem.max())
#		plt_axes('x', 'both', [10,1])
#		plt.ylabel('[A.U.]', size=16)
#		plt.ylim(-vmax, vmax)
#		plt.title('Imaginary part', size=16)
#		plt.suptitle('Fitted fringe of baseline='+strbl+' @ '+freqstr+'MHz by FitVis()', fontsize=20)
#		plt.savefig(self.outdir+'vis_FitVis_'+strbl+'_'+freqstr+'MHz_'+str(self.nsigma)+'.png')
#		plt.close()
#		return [self.Ampv, self.Sigmav, self.Phaseaddv]
#
#
#
#	def FitExp( self, Nprocess=None ):
#		if ('Timeb' not in self.__dict__.keys()) : self.FitBeam()
#		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
#		freq = self.antarray.Ant.freq
#		inttime = self.antarray.Ant.inttime
#		x0 = np.arange(self.vis.shape[0])*inttime  # second
#		#--------------------------------------------------
#		# Move time axis to -1
#		vis = ArrayAxis(self.vis.data, 0, -1, 'move')
#		shape = vis.shape  # (512,528,3600)
#		# Reshape to 2D
#		vis = vis.reshape(np.prod(shape[:2]), shape[-1])
#		#--------------------------------------------------
#		pool = jp.PoolFor(0, len(vis), Nprocess)
#		send = (vis, self.Timeb.flatten())
#		bcast = (x0, self.Dec, freq, bl, self.showprogress, shape[1], self.antarray.Ant.lonlat[1])
#		pf = pool.map_async(_DoMultiprocess_FitExp, send, bcast)
#		#--------------------------------------------------
#	#	self.Phaseadde = np.concatenate(pf) %(2*np.pi)
#	#	self.Phaseadde = self.Phaseadde.reshape(shape[:2])
#		pf = np.concatenate(pf).reshape(shape[:2])
#		self.Phaseadde = pf %(2*np.pi)
#		del pool, bcast, vis
#		#--------------------------------------------------
#		np.save(self.outdir+'Phaseadde_'+str(self.nsigma)+'.npy', self.Phaseadde)
#		#--------------------------------------------------
#		# Plot
#		if (plotnv is None) : 
#			bl = abs(bl[:,0]) # longest East-West baseline
#			nv = np.where(bl==bl.max())[0][0]
#		else : nv = plotnv
#		nvbl = np.arange(len(self.antarray.Blorder.blorder))[self.antarray.visorder][nv]
#		bl = self.antarray.Blorder.Order2Bl(nvbl)
#		strbl = str(bl[0])+'-'+str(bl[1])
#
#		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
#		if (plotfreq is None) : plotfreq = freq.mean()
#		freqstr = str(int(round(plotfreq)))
#		nf = abs(freq-plotfreq)
#		nf = np.where(nf==nf.min())[0][0]
#		vis = self.vis[:,nf,nv]  # real+imag+masked
#		vmax = abs(vis).max()*1.05
#		A = abs(vis.data)
#		#--------------------
#		def func(x, p) : 
#			beta = jp.Circ2Sph((x-self.Timeb[nf,nv])/3600.*15*np.pi/180, self.Dec*np.pi/180)
#			angle = 2*np.pi/(300./freq[nf]) * ( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((self.Dec-self.antarray.Ant.lonlat[1])*np.pi/180) ) + p[0]
#			return A*np.exp(1j*angle)  # complex
#		#--------------------
#		visf = func(x0, [self.Phaseadde[nf,nv]])
#		plt.figure(figsize=(17,6))
#		plt.subplot(1,2,1)
#		plt.plot(self.timem, vis.real, 'b-', label='data')
#		plt.plot(self.timem, visf.real, 'r-', label='Fitting')
#		plt.legend(fontsize=16)
#		plt.xlabel('$t$ [min]', size=16)
#		plt.xlim(self.timem.min(), self.timem.max())
#		plt_axes('x', 'both', [10,1])
#		plt.ylabel('[A.U.]', size=16)
#		plt.ylim(-vmax, vmax)
#		plt.title('Real part', size=16)
#		plt.subplot(1,2,2)
#		plt.plot(self.timem, vis.imag, 'b-', label='data')
#		plt.plot(self.timem, visf.imag, 'r-', label='Fitting')
#		plt.legend(fontsize=16)
#		plt.xlabel('$t$ [min]', size=16)
#		plt.xlim(self.timem.min(), self.timem.max())
#		plt_axes('x', 'both', [10,1])
#		plt.ylabel('[A.U.]', size=16)
#		plt.ylim(-vmax, vmax)
#		plt.title('Imaginary part', size=16)
#		plt.suptitle('Fitted fringe of baseline='+strbl+' @ '+freqstr+'MHz by FitExp()', fontsize=20)
#		plt.savefig(self.outdir+'vis_FitExp_'+strbl+'_'+freqstr+'MHz_'+str(self.nsigma)+'.png')
#		plt.close()



	def Plot( self, dyDeff=None, dyLew=None, Nprocess=None ) : 
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliPhase.Plot: start @', starttime
			if (self.vistype == 'auto') : 
				print '    CaliPhase.vistype='+self.vistype
				endtime = jp.Time(1)
				costtime = jp.Time(starttime, endtime)
				tottime = jp.Time(self.starttime, endtime)
				print 'CaliPhase.Plot:  end  @', endtime
				print 'Cost:', costtime+'     Total:', tottime+'\n'
		outdir = jp.Outdir((None,'file'), (0,'file'), (0,'func'))
		#--------------------------------------------------
		nbl = np.arange(len(self.antarray.visorder))
		bl = self.antarray.Blorder.blorder[self.antarray.visorder]
		strbl = []
		for i in xrange(len(nbl)) : 
			strbl.append(str(bl[i,0])+'-'+str(bl[i,1]))
		#--------------------------------------------------
		self._PlotAmp( outdir, strbl, Nprocess)
		self._PlotDeff(outdir, strbl, Nprocess, dyDeff)
		self._PlotLew( outdir, strbl, Nprocess, dyLew)
		self._PlotPhaseadd(outdir, strbl, Nprocess)
		if (self.verbose) : 
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliPhase.Plot:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'



	def _PlotAmp( self, outdir, strbl, Nprocess ) : 
		freq = self.antarray.Ant.freq
		freqminmax = [freq.min(), freq.max()]
		freq = freq[self.antarray.freqorder]
		func, amp, color = ['FitBeam', 'FitVis'], [], ['b','r']
		try    : amp.append(self.Ampb.T)  # (Nv,Nf)
		except : amp.append(None)
		try    : amp.append(self.Ampv.T)
		except : amp.append(None)
		if (amp[0] is None and amp[1] is None) : 
			if (self.verbose) : print '    NOT exist  self.Ampb, self.Ampv'
			return
		if (self.verbose) : print '    Plotting    Amp     ......'
		send = tuple(amp) + (strbl,)
		bcast = (freq, func, color, outdir, 'amp', self.nsigma, freqminmax)
		if (Nprocess <= 1) : 
			iterable = [[(0,len(strbl)),None], send, bcast]
			_DoMultiprocess_Plot( iterable )
		else : 
			pool = jp.PoolFor(0, len(strbl), Nprocess)
			pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotDeff( self, outdir, strbl, Nprocess, dyDeff ) : 
		freq = self.antarray.Ant.freq
		freqminmax = [freq.min(), freq.max()]
		freq = freq[self.antarray.freqorder]
		func, Deff, color = ['FitBeam', 'FitVis'], [], ['b','r']
		try    : Deff.append(300/(self.Sigmab.T*self.fwhm2sigmafactor*freq[None,:]))
		except : Deff.append(None)
		try    : Deff.append(300/(self.Sigmav.T*self.fwhm2sigmafactor*freq[None,:]))
		except : Deff.append(None)
		if (Deff[0] is None and Deff[1] is None and self.verbose) : 
			if (self.verbose) : print '    NOT exist  self.Sigmab, self.Sigmav'
			return
		if (self.verbose) : print '    Plotting    Deff    ......'
		send = tuple(Deff) + (strbl,)
		antform = jp.npfmt(self.antarray.Ant.antform)[0]
		bcast = (freq, func, color, outdir, 'deff', self.nsigma, freqminmax, dyDeff, self.fwhm2sigmafactor, antform)
		if (Nprocess <= 1) : 
			iterable = [[(0,len(strbl)), None], send, bcast]
			_DoMultiprocess_Plot( iterable )
		else : 
			pool = jp.PoolFor(0, len(strbl), Nprocess)
			pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotLew( self, outdir, strbl, Nprocess, dyLew ) : 
		freq = self.antarray.Ant.freq
		freqminmax = [freq.min(), freq.max()]
		freq = freq[self.antarray.freqorder]
		func, Lew, color = ['FitPhase-m0','FitPhase-m1','FitPhase-m2','FitVis','FitExp'], [], ['b','r','c','k','y']
		try    : Lew.append(self.Lewp0.T)  # (Nv,Nf)
		except : Lew.append(None)
		try    : Lew.append(self.Lewp1.T)
		except : Lew.append(None)
		try    : Lew.append(self.Lewp2.T)
		except : Lew.append(None)
		if (Lew[0] is None and Lew[1] is None) : 
			if (self.verbose) : print '    NOT exist  self.Lewp, self.Lewv'
			return
		if (self.verbose) : print '    Plotting    Lew     ......'
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		send = tuple(Lew) + (strbl,)
		bcast = (freq, func, color, outdir, 'lew', self.nsigma, freqminmax, self.fitLew, bl, dyLew)
		if (Nprocess <= 1) : 
			iterable = [[(0,len(strbl)), None], send, bcast]
			_DoMultiprocess_Plot( iterable )
		else : 
			pool = jp.PoolFor(0, len(strbl), Nprocess)
			pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotPhaseadd( self, outdir, strbl, Nprocess ) : 
		freq = self.antarray.Ant.freq
		freqminmax = [freq.min(), freq.max()]
		freq = freq[self.antarray.freqorder]
		func, phaseadd, color = ['FitPhase-m0','FitPhase-m1','FitPhase-m2','FitVis','FitExp'], [], ['b','r','c','k','y']
		try    : phaseadd.append(self.Phaseaddp0.T*180/np.pi) # deg
		except : phaseadd.append(None)
		try    : phaseadd.append(self.Phaseaddp1.T*180/np.pi) # deg
		except : phaseadd.append(None)
		try    : phaseadd.append(self.Phaseaddp2.T*180/np.pi) # deg
		except : phaseadd.append(None)
		try    : phaseadd.append(self.Phaseaddv.T*180/np.pi)
		except : phaseadd.append(None)
		try    : phaseadd.append(self.Phaseadde.T*180/np.pi)
		except : phaseadd.append(None)
		if (phaseadd[0] is None and phaseadd[1] is None and phaseadd[2] is None) : 
			if (self.verbose) : print '    NOT exist  self.Phaseaddp, self.Phaseaddv, self.Phaseadde'
			return
		if (self.verbose) : print '    Plotting  Phaseadd  ......'
		bl =self.antarray.Blorder.baseline[self.antarray.visorder]
		send = tuple(phaseadd) + (strbl,)
		bcast = (freq, func, color, outdir, 'phaseadd', self.nsigma, freqminmax, self.fitLew, bl)
		if (Nprocess <= 1) : 
			iterable = [[(0,len(strbl)), None], send, bcast]
			_DoMultiprocess_Plot( iterable )
		else : 
			pool = jp.PoolFor(0, len(strbl), Nprocess)
			pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def PlotFringe( self, nf, nv, array, label, color, lw=None ):
		'''
		Example:
			CaliPhase.PlotFringe(0, np.arange(antarray.visorder.size), caliphase.vis, ['Real','Imag','Abs'], ['b','r','c'], None)

		nf, nv:
			Used to array[:,nf,nv], can be int or np.array([,int])
			Which to be plotted
			nf==None: center frequency
			nv=None: longest East-West baseline

		array: 
			vis to be plotted
			(1) ndarray
			(2) [ndarray, ndarray, ...], in this case, label, color, lw should has the same size

		label: 
			plt.plot(x, y, label=label)
			label must be str pair [label1, label2], 1=>real, 2=>imag
			For array=[array1, array2], label=[[label11,label12], [label21,label22]]

		color:
			color must be str pair [color1, color2], 1=>real, 2=>imag

		lw:
			linewidth, must be int, one int for both real+imag
			can be lw=None
		'''
		outdir = jp.Outdir((None,'file'), (0,'file'), (0,'func'))
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : 
			print 'CaliPhase.PlotFringe'
		if (self.timem[-1]-self.timem[0] <= 60) : xmajor = 5
		else : xmajor = 10
		#--------------------------------------------------

		istype = jp.IsType()
		if (not (istype.islist(array) or istype.istuple(array))) : 
			if (lw is None) : lw = 1 
			array,label,color,lw = [array],[label],[color],[lw]
		elif (lw is None) : lw = [1,1] 
		if (not len(array)==len(label)==len(color)==len(lw)) : 
			jp.Raise(Exception, 'Error: len(array)='+str(len(array))+', len(label)='+str(len(label))+', len(color)='+str(len(color))+', len(lw)='+str(len(lw))+' are NOT equal')
		#--------------------------------------------------

		# Plot
		nf, strfreq, nv, strblo, strbli, blname = self.antarray.Plotnfnv(nf, nv)
		for j in xrange(nf.size) : 
			for k in xrange(nv.size) : 

				vmax, vmin = [], []
				if (self.vistype == 'cross') : plt.figure(figsize=(25,6))
				for i in xrange(len(array)) : 

					vmax.append( abs(array[i][:,nf[j],nv[k]]).max()*1.05 )
					vmin.append( abs(array[i][:,nf[j],nv[k]]).min()*0.95 )
					scale = jp.SciNot(vmax[-1])[1]-2
					vmax[-1] /= 10.**scale
					vmin[-1] /= 10.**scale
					#------------------------------
					if (self.vistype == 'cross') : 
						plt.subplot(1,3,1)
						plt.plot(self.timem, array[i].real[:,nf[j],nv[k]]/10.**scale, color=color[i][0], lw=lw[i], label=label[i][0])
					else : 
						plt.plot(self.timem, array[i][:,nf[j],nv[k]]/10.**scale, color=color[i][0], lw=lw[i], label=label[i][0])
					#------------------------------
					if (i == len(array)-1) : 
						plt.legend()
						plt.xlabel('time [min]', size=16)
						plt.xlim(self.timem.min(), self.timem.max())
						plt_axes('x', 'both', [xmajor,1])
						plt.ylabel('[A.U.]', size=16)
						vmax = np.array(vmax).max()
						vmin = np.array(vmin).min()
						if (self.vistype == 'cross') : plt.ylim(-vmax, vmax)
						else : plt.ylim(vmin, vmax)
					if (self.vistype != 'cross') : continue
					#------------------------------
					plt.subplot(1,3,2)
					plt.plot(self.timem, array[i].imag[:,nf[j],nv[k]]/10.**scale, color=color[i][1], lw=lw[i], label=label[i][1])
					if (i == len(array)-1) : 
						plt.legend()
						plt.xlabel('time [min]', size=16)
						plt.xlim(self.timem.min(), self.timem.max())
						plt_axes('x', 'both', [xmajor,1])
						plt.ylabel('[A.U.]', size=16)
						plt.ylim(-vmax, vmax)
				#--------------------------------------------------
					plt.subplot(1,3,3)
					plt.plot(self.timem, abs(array[i][:,nf[j],nv[k]]/10.**scale), color=color[i][2], lw=lw[i], label=label[i][2])
					if (i == len(array)-1) : 
						plt.legend()
						plt.xlabel('time [min]', size=16)
						plt.xlim(self.timem.min(), self.timem.max())
						plt_axes('x', 'both', [xmajor,1])
						plt.ylabel('[A.U.]', size=16)
						plt.ylim(0, vmax)
				#--------------------------------------------------
				if ( jp.SysFrame(0,2)[-1][-2]=='' ) :  
					if (self.vistype=='cross' and os.path.exists(outdir+'vis_'+strblo[k]+'_'+strfreq[j]+'MHz_'+str(self.nsigma)+'.png')) : 
						if (self.verbose) : print '    Exist cross '+strblo[k]+'='+strbli[k]+' @ '+strfreq[j]+'MHz'
						plt.close()
						continue
					elif (self.vistype=='auto' and outdir+'vis_'+strbl+'_'+strfreq[j]+'MHz_'+str(self.nsigma)+'.png') : 
						if (self.verbose) : print '    Exist auto '+strbl+' @ '+strfreq[j]+'MHz'
						plt.close()
						continue
				#--------------------------------------------------
				if (self.vistype == 'cross') : 
					plt.suptitle('Fringe of baseline='+strblo[k]+'='+strbli[k]+' @ '+strfreq[j]+'MHz', fontsize=20)
					plt.savefig(outdir+'vis_'+strblo[k]+'_'+strfreq[j]+'MHz_'+str(self.nsigma)+'.png')
					if (self.verbose) : print '    Plotting cross '+strblo[k]+strbli[k]+' @ '+strfreq[j]+'MHz'
				else : 
					strbl = strblo[k][:strblo[k].find('-')]
					plt.title('Auto-correlation of channel='+strbl+' @ '+strfreq[j]+'MHz', fontsize=16)
					plt.savefig(outdir+'vis_'+strbl+'_'+strfreq[j]+'MHz_'+str(self.nsigma)+'.png')
					if (self.verbose) : print '    Plotting auto '+strbl+' @ '+strfreq[j]+'MHz'
				plt.close()
				if (self.vistype != 'cross') : return
				#--------------------------------------------------

				# Plot phase
				label = ['Phase', 'Smoothed']
				markersize = [4, 3]
				colorls = ['bo', 'ro']
				for i in xrange(len(array)) : 
					phase = np.angle(array[i][:,nf[j],nv[k]])
					phase = (phase*180/np.pi) %360
					plt.plot(self.timem, phase, colorls[i], markersize=markersize[i], label=label[i])
				plt.legend(fontsize=10)
				plt.xlabel('time [min]', size=16)
				plt.xlim(self.timem.min(), self.timem.max())
				plt_axes('x', 'both', [xmajor,1])
				plt.ylim(0, 360)
				plt.ylabel('[degree]', size=16)
				plt.title('Total phase of baseline='+strblo[k]+strbli[k]+' @ '+strfreq[j]+'MHz', size=13)
				plt.savefig(outdir+'phase_'+strblo[k]+'_'+strfreq[j]+'MHz_'+str(self.nsigma)+'.png')
				plt.close()
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : print







