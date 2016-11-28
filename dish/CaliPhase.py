import os
import time
from copy import deepcopy as dcopy
from AntArray import *
from Masking import *
import jizhipy as jp
from jizhipy.Plot import *

#import scipy.signal as spsn
#from Sph2Circ import *
#from Plot import *
#from PoolFor import *
#from ResetMasked import *
#from Smooth import *
#from FuncFit import *

##################################################



def _DoMultiprocess_FitBeam( iterable ) : 
	n1, n2  = iterable[0]
	beam = iterable[1]
	timem, s0, Dec, verbose = iterable[2]
	if (verbose) : progressbar = jp.ProgressBar('    pid='+str(os.getpid())+':', len(beam), False)
	#--------------------------------------------------
	def func(x, p) :  # x is time, not angle
		beta = jp.Circ2Sph((x-p[1])/60.*15*np.pi/180, Dec*np.pi/180)
		y = p[0] * np.exp(-beta**2 /2 /p[2]**2)
		return y
	#--------------------------------------------------
	pf = []
	for i in xrange(n2-n1) : 
		if (verbose) : progressbar.Progress()
		p0 = [beam[i].max(), timem.mean(), s0]
		pf.append( jp.FuncFit(func, timem, beam[i], p0, warning=False)[0] )
	return jp.npfmt(pf)



def _DoMultiprocess_FitPhase( iterable ) : 
	n1, n2 = iterable[0]
	phase, t0list = iterable[1]
	timem, Dec, freq, bl, verbose, Nf, lat, fitLew = iterable[2]
	if (verbose) : progressbar = jp.ProgressBar('    pid='+str(os.getpid())+':', len(phase), False)
	pf = []
	#--------------------------------------------------
	for i in xrange(n2-n1) : 
		t0 = t0list[i] #@
		if (verbose) : progressbar.Progress()
		nv, nf = (i+n1)/Nf, (i+n1)%Nf
		#--------------------------------------------------
		def funcNO(x, p) : 
			beta =jp.Circ2Sph((x-t0)/60.*15*np.pi/180, Dec*np.pi/180)
			angle = 2*np.pi/(300./freq[nf]) *( bl[nv,0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
			return np.sin(angle)
		pf.append( list(jp.FuncFit(funcNO, timem, np.sin(phase[i]), [1])[0]) )
		#--------------------------------------------------
		if (fitLew == 1) : 
			def funcYES(x, p) : 
				beta =jp.Circ2Sph((x-t0)/60.*15*np.pi/180, Dec*np.pi/180)
				angle = 2*np.pi/(300./freq[nf]) *( p[0]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + pf[-1][0]
				return np.sin(angle)
			pf[-1] += list(jp.FuncFit(funcYES, timem, np.sin(phase[i]), bl[nv,0])[0])
		elif (fitLew == 2) : 
			def funcBOTH(x, p) : 
				beta = jp.Circ2Sph((x-t0)/60.*15*np.pi/180, Dec*np.pi/180)
				angle = 2*np.pi/(300./freq[nf]) *( p[1]*np.sin(beta) - bl[nv,1]*np.cos(beta)*np.sin((Dec-lat)*np.pi/180) ) + p[0]
				return np.sin(angle)
			pf[-1] = jp.FuncFit(funcBOTH, timem, np.sin(phase[i]), [pf[-1][0], bl[nv,0]])[0]
		#--------------------------------------------------
		elif (fitLew == 0) : pf[-1] += [bl[nv,0]]
	return jp.npfmt(pf)



def _DoMultiprocess_FitVis( iterable ) : 
	n1, n2 = iterable[0]
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
	n1, n2 = iterable[0]
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
	n1, n2 = iterable[0]
	array = iterable[1][:-1]
	strbl = iterable[1][-1]
	freq, func, color, outdir, which, nsigma = iterable[2][:6]
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
		plt.xlim(int(round(freq.min())), int(round(freq.max())))
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


#	def __init__( self, antarray=None, masking=None, caligain=None, nsigma=4, fwhm2sigmafactor=2.287, fitLew=2, plotnv=None, plotfreq=None, Nprocess=None, verbose=True, outdir=None ) : 
	def __init__( self, antarray=None, masking=None, caligain=None, fwhm2sigmafactor=2.287, Nprocess=None, verbose=True, outdir=None ) : 
		'''
		antarray: instance of class:AntArray
		masking:  instance of class:Masking
		caligain: instance of class:CaliGain

		fwhm2sigmafactor:
			sigma = lambda /Deff /fwhm2sigmafactor

#		nsigma: 
#			times of the sigma of the beam/fringe
#
#		fitLew:
#			=0, 1, 2
#			Different metohds to fit Lew
#
#		plotfreq:
#			in MHz
#
#		plotnv:
#			AntArray.visorder[plotnv]
#			plotnv is the number/index/order of visorder, NOT the number/index/order of AntArray.Blorder.blorder
#
#		If ==None (except Nprocess), don't do __init__() it
		'''
		self.starttime = jp.Time(1)
		class _Params( object ) : pass
		self.Params = _Params()
		self.Params.init = {'nsigma':nsigma, 'fwhm2sigmafactor':fwhm2sigmafactor, 'fitLew':fitLew, 'Nprocess':Nprocess, 'verbose':verbose}
		if (verbose is not None) : self.verbose = verbose
		if (self.verbose) : print '\n'
		self.Nprocess = jp.NprocessCPU(Nprocess, verbose)[0]
		self.plotfreq, self.plotnv = plotfreq, plotnv
		#--------------------------------------------------
		if (nsigma is not None) : 
			try : self.nsigma = float(nsigma)
			except : self.nsigma = 4
			if (self.nsigma == int(self.nsigma)) : self.nsigma = int(self.nsigma)
		#--------------------------------------------------
		if (fwhm2sigmafactor is not None) : 
			try : self.fwhm2sigmafactor = float(fwhm2sigmafactor)
			except : self.fwhm2sigmafactor = 2.287
		#--------------------------------------------------
		if (fitLew is not None) : 
			try : self.fitLew = int(fitLew)
			except : self.fitLew = 2
		#--------------------------------------------------
#		if (antarray is not None) : 
#			if (type(antarray) not in [tuple, list]) : 
#				antarray = [antarray]
#			count = np.zeros(len(antarray), int)
#			nhdf5 = []
#			for i in xrange(len(antarray)) : 
#				if (antarray[i].vistype not in ['cross1','cross2']) :
#					count[i] = 1
#				nhdf5.append(antarray[i].Hdf5.nhdf5)
#			nhdf5transit = antarray[0].Hdf5.nhdf5transit[-1]
#			if (nhdf5transit not in nhdf5) : jp.Raise(Exception, 'antarray.Hdf5.nhdf5='+str(nhdf5)+' is NOT the files which contains the transit source (nhdf5='+str(nhdf5transit)+')')
#			self.antarray = antarray  # share address, not copy
#			if (count.sum() == 0) : self.vistype = 'cross'
#			else : self.vistype = 'auto'
#		#--------------------------------------------------
#		#	self.antarray.SelectVisType()  # Ensure cross1
#		if (masking is not None) : 
#			if (type(masking) not in [tuple, list]) : 
#				masking = [masking]
#			self.masking = masking
#		#--------------------------------------------------
#		if (caligain is not None) : 
#			if (type(caligain) not in [tuple, list]) : 
#				caligain = [caligain]
#			self.caligain = caligain
#		#--------------------------------------------------
		if (antarray is not None) : 
			self.antarray = antarray
			try : a, b = antarray.Hdf5.transittime, antarray.Hdf5.nhdf5transit
			except : jp.Raise(Exception, "NOT exist fo['transitsource'], you can set it by hand with AntArray.Transitsource()")
		#--------------------------------------------------
		if (masking is not None) : self.masking = masking
		if (caligain is not None) : self.caligain = caligain
		self.outdir = jp.Outdir((None,'file'), (0,'file'))
		if (outdir is not None) : self.outdir = jp.Mkdir(self.outdir+outdir)



	def RADec( self, RAORsourcename=None, Dec=None ) : 
		'''
		Dec: Dec of the data/antenna pointing/calibration source
		RA : RA  of the calibration source 
		angle in degree
		'''
		istype = jp.IsType()
		if (istype.isstr(RAORsourcename)) : 
			brightsource = jp.BrightSource()
			self.RA, self.Dec = brightsource.RADec(RAORsourcename)
		else : 
			if (RA  is not None) : self.RA  = RA
			if (Dec is not None) : self.Dec = Dec
		self.Params.RADec = {'RA':self.RA, 'Dec':self.Dec}



	def Fringe( self, whichtransitsource=-1 ) : 
		'''
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
#		self.yshift = []
		inttime = self.antarray.Ant.inttime
		freq = self.antarray.Ant.freq
		antform = jp.npfmt(self.antarray.Ant.antform)[0]
#		Nt = self.antarray[0].vis.shape[0]  # each hdf5 must have the same number of points
#		timerange1, timerange2, ncount = -10, Nt+10, 0
		# Here we assume: fringe must be in 1 or 2 hdf5 files, not more than 2
		#--------------------------------------------------

		pixsource = int(self.antarray.Hdf5.transittime / inttime) - self.antarray.Hdf5.nhdf5[0] * self.antarray.Ant.N0
		#--------------------------------------------------

		sigma = 300/freq.mean()/antform/self.fwhm2sigmafactor
		sigma = jp.Sph2Circ(sigma, self.Dec*np.pi/180)  # rad
		# Total pixels of the fringe in 1 sigma
		pix1 = int(round(24*3600/(2*np.pi)*sigma/inttime))

		pixpair10 = [pixsource-10*pix1/2, pixsource+10*pix1/2]
		if (pixpair10[0] < 0) : pixpair10[0] = 0
		if (pixpair10[1] > self.antarray.vis.shape[0]) : pixpair10[1] = self.antarray.vis.shape[0]
#		timem10 = self.antarray.timem[pixpair10[0], pixpair10[1]]

		if (self.vistype == 'cross') : self._PlotFringe(visTEN.data, ['Real part', 'Imaginary part'], ['b', 'r'])
		



















		timerangeONE = int(round(24*3600/(2*np.pi)*sigma/inttime))
		# 10 sigma range

		nsigmaTEN = 15
		timerangeTEN1, timerangeTEN2 = transittime-nsigmaTEN*timerangeONE/2, transittime+nsigmaTEN*timerangeONE/2+1  # pixels
		# self.nsigma
		timerangeSELF1, timerangeSELF2 = transittime-self.nsigma*timerangeONE/2, transittime+self.nsigma*timerangeONE/2+1  # pixels

		trTEN  = np.arange(timerangeTEN1, timerangeTEN2)
		trSELF = np.arange(timerangeSELF1, timerangeSELF2)

		trTENfile, trTENlocal = trTEN / self.antarray[0].Nt, trTEN % self.antarray[0].Nt
		trSELFfile, trSELFlocal = trSELF / self.antarray[0].Nt, trSELF % self.antarray[0].Nt

		nhdf5TEN = np.arange(trTENfile[0], trTENfile[-1]+1)
		nhdf5SELF = np.arange(trSELFfile[0], trSELFfile[-1]+1)

		timerangeTEN, timerangeSELF = [], []
		for i in xrange(nhdf5TEN.size) : 
			ntmp = trTENlocal[trTENfile==nhdf5TEN[i]]
			timerangeTEN.append([ntmp[0], ntmp[-1]+1])
		for i in xrange(nhdf5SELF.size) : 
			ntmp = trSELFlocal[trSELFfile==nhdf5SELF[i]]
			timerangeSELF.append([ntmp[0], ntmp[-1]+1])
		timerangeTEN, timerangeSELF = np.array(timerangeTEN), np.array(timerangeSELF)
		#--------------------------------------------------

		# Plot TEN
		visTEN, maskTEN = [], []
		antarray = AntArray()
		antarray.__dict__ = self.antarray[0].__dict__
		antarray.verbose = False
		for i in xrange(nhdf5TEN.size) : 
			antarray.WhichHdf5(nhdf5TEN[i])
			masking = Masking(antarray, verbose=False)
			try : masking.MaskNoiseSource(self.masking[0].noisesource.pixstart, self.masking[0].noisesource.pixlength, self.masking[0].noisesource.pixperiod)
			except : pass
			vis=antarray.vis[:,antarray.freqorder,antarray.visorder]
			if (antarray.freqorder.size == 1) : vis = vis[:,None]
			if (antarray.visorder.size == 1) : vis = vis[:,:,None]
			vis[masking.mask] = masking.maskvalue
			visTEN.append(vis[timerangeTEN[i,0]:timerangeTEN[i,1]])
			maskTEN.append(masking.mask[timerangeTEN[i,0]:timerangeTEN[i,1]])
		visTEN = np.concatenate(visTEN, 0)
		maskTEN = np.concatenate(maskTEN, 0)
		visTEN = np.ma.MaskedArray(visTEN, maskTEN)
		maskTEN = antarray = masking = vis = 0 #@
		self.yshift10 = visTEN.mean(0)  #@#@
		visTEN -= self.yshift10

		self.timem = trTEN * self.antarray[0].Ant.inttime /60
		if (self.vistype == 'cross') : self._PlotFringe(visTEN.data, ['Real part', 'Imaginary part'], ['b', 'r'])
		else : self._PlotFringe(visTEN.data, ['Auto'], ['b'])
		#--------------------------------------------------

		np.save(self.outdir+'timem_CygA_1', self.timem)
		print visTEN.data.shape
		np.save(self.outdir+'vis_CygA_1', visTEN.data.flatten())

		jp.Raise()










#			visi = self.antarray[i].vis[:,self.antarray[i].freqorder,self.antarray[i].visorder]
#			nf, nv = self.antarray[i].freqorder.size, self.antarray[i].visorder.size
#			if (nf == 1) : visi = visi[:,None]
#			if (nv == 1) : visi = visi[:,:,None]
#			#--------------------------------------------------
#
#			try : 
#				visi[self.masking[i].mask]= self.masking[i].maskvalue
#				mask.append(self.masking[i].mask)
#			except : 
#				lx = (timerangeTEN[i,1]-timerangeTEN[i,0],)
#				mask.append( np.zeros(lx+visi.shape[1:], bool) )
#			#--------------------------------------------------
#
#			try : 
#				x = jp.Edge2Center(self.caligain[i].nsplit)
#				y = self.caligain[i].gaintmean
#				xnew = np.arange(visi.shape[0])
#				gaintmean = []
#				for j in xrange(y.shape[-1]) : 
#					gaintmean.append(jp.Interp1d(x, y[:,0,j], new))
#				visi /= np.array(gaintmean)[:,None,:]
#			except : pass
#			#--------------------------------------------------
#
#			self.yshift.append( visi.mean(0) )
#			self.vis.append( visi[timerangeTEN[i,0]:timerangeTEN[i,1]] )
#		#--------------------------------------------------

		if (len(self.vis) == 1) : 
			self.vis = self.vis[0]
			mask = mask[0]
		else : 
			self.vis = np.concatenate(self.vis, 0)
			mask = np.concatenate(mask, 0)
		#--------------------------------------------------

		if (self.vistype == 'cross') : 
			if (len(self.yshift) == 2) : self.yshift = ((self.yshift[0]+self.yshift[1])/2.)[None,:,:]
			else : self.yshift = self.yshift[0][None,:,:]
			self.vis -= self.yshift
		else : self.vis = self.vis.real
		#--------------------------------------------------

		print self.vis.shape, mask.shape
		jp.Raise()
		self.vis   = np.ma.MaskedArray(self.vis, mask)  #@#@
		Npix = (timerange[:,1]-timerange[:,0]).sum()
		self.timem = np.arange(timerange[0,0], timerange[0,0]+Npix) *inttime/60. #@#@ min
		mask = 0 #@
		#--------------------- Plot -----------------------
		if (self.vistype == 'cross') : self._PlotFringe(self.vis, ['Real part', 'Imaginary part'], ['b', 'r'])
		else : self._PlotFringe(self.vis, ['Auto'], ['b'])
		#--------------------------------------------------
		if (nsigma < self.nsigma) : 
			if (nsigma == int(nsigma)) : nsigma = int(nsigma)
			dn = int(1.*nsigma/self.nsigma * self.timem.size)
			n1, n2 = self.timem.size/2-dn/2, self.timem.size/2+dn/2
			self.timemlarger = self.timem.copy()
			self.vislarger   = self.vis.copy()
			self.timem = self.timem[n1:n2]
			self.vis   = self.vis[n1:n2]
			np.save(self.outdir+'timemlarger_minute_'+str(self.nsigma), self.timemlarger)
			np.save(self.outdir+'timem_minute_'+str(nsigma), self.timem)
			self.nsigma = nsigma
			if (self.vistype == 'cross') : self._PlotFringe(self.vis, ['Real part', 'Imaginary part'], ['b', 'r'])
			else : self._PlotFringe(self.vis, ['Auto'], ['b'])
		else : 
			self.timemlarger = self.timem
			self.vislarger   = self.vis
			np.save(self.outdir+'timem_minute_'+str(self.nsigma), self.timem)
			np.save(self.outdir+'timemlarger_minute_'+str(self.nsigma), self.timemlarger)
		#!!!!! vis,       timem      : used to fit  !!!!!
		#!!!!! vislarger, timemlarger: used to plot !!!!!
		#--------------------------------------------------
		if (self.verbose) : 
			print '    timem_minute_'+str(nsigma)+'.npy, timemlarger_minute_'+str(self.nsigma)+'.npy  saved'
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliPhase.Fringe:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'




	def Smooth( self, timetimes=100 ) : 
		''' Smooth(self.vis, 0, 3, timetimes) '''
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliPhase.Smooth: start @', starttime
		array = [self.vis]
		mask = self.vis.mask.copy()
		self.vis = jp.Smooth(self.vis.data, 0, 3, timetimes, Nprocess=self.Nprocess)
		self.vis = np.ma.MaskedArray(self.vis, mask)
		array.append(self.vis)
		if (self.vistype == 'cross') : self._PlotFringe(array, [['Real part','Imaginary part'],['Smoothed','Smoothed']], [['b','r'],['m','g']], [1,3])
		else : self._PlotFringe(array, [['Auto-corr'],['Smoothed']], [['b'],['r']], [1,3])
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
				print '    CaliPhase.vistype='+self.vistype
				endtime = jp.Time(1)
				costtime = jp.Time(starttime, endtime)
				tottime = jp.Time(self.starttime, endtime)
				print 'CaliPhase.FitBeam:  end  @', endtime
				print 'Cost:', costtime+'     Total:', tottime+'\n'
		# Initial guess
		inttime = self.antarray[0].Ant.inttime
		antform = jp.npfmt(self.antarray[0].Ant.antform)[0]
		s0 = 300/self.antarray[0].Ant.freq.mean()/antform/2.287
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
		pf = np.concatenate(pf)  # (nv*nf, 3)
		pf = pf.reshape(shape[:2]+(pf.shape[-1],)).T  # (3,nf,nv)
		self.Ampb, self.Timeb, self.Sigmab = pf
		del pool, beam
		self.Ampb = abs(self.Ampb)
		self.Sigmab = abs(self.Sigmab)
		#--------------------------------------------------
		np.save(self.outdir+'Ampb_'+str(self.nsigma)+'.npy', self.Ampb)
		np.save(self.outdir+'Timeb_'+str(self.nsigma)+'.npy',self.Timeb)
		np.save(self.outdir+'Sigmab_'+str(self.nsigma)+'.npy', self.Sigmab)
		bl =self.antarray[0].Blorder.baseline[self.antarray[0].visorder]
		self.Phasens = 2*np.pi/300*self.antarray[0].Ant.freq[:,None]*(-bl[:,1][None,:])*np.sin((self.Dec-self.antarray[0].Ant.lonlat[1])*np.pi/180)
		np.save(self.outdir+'Phasens.npy', self.Phasens)
		# vis * exp(-1j*Phasens) * exp(-1j*Phaseadd)
		if (self.verbose) : 
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliPhase.FitBeam:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'



	def FitPhase( self, fitLew=None ) : 
		'''
		(1) Because there is dPhase, fringe is not sensitive enough to Lew, we need to first fit the initial guess of dPhase.
		(2) Also because fringe is not sensitive enough to Lew, different range of fringe (nsigma=1, 2, 3, ...) have very different Lew.
		'''
		if (fitLew is not None) : 
			try : self.fitLew = int(fitLew)
			except : self.fitLew = 2
		#--------------------------------------------------
		if ('Timeb' not in self.__dict__.keys()) : self.FitBeam()
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliPhase.FitPhase: start @', starttime
			if (self.vistype == 'auto') : 
				print '    CaliPhase.vistype='+self.vistype
				endtime = jp.Time(1)
				costtime = jp.Time(starttime, endtime)
				tottime = jp.Time(self.starttime, endtime)
				print 'CaliPhase.FitPhase:  end  @', endtime
				print 'Cost:', costtime+'     Total:', tottime+'\n'
		bl =self.antarray[0].Blorder.baseline[self.antarray[0].visorder]
		freq = self.antarray[0].Ant.freq
		inttime = self.antarray[0].Ant.inttime
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
		bcast = (self.timem, self.Dec, freq, bl, self.verbose, shape[1], self.antarray[0].Ant.lonlat[1], self.fitLew)
		pf = pool.map_async(_DoMultiprocess_FitPhase, send, bcast)
		print
		#--------------------------------------------------
		pf = np.concatenate(pf)  # (nv*nf, 2)
		pf = pf.reshape(shape[:2]+(pf.shape[-1],)).T  # (2,nf,nv)
		self.__dict__['Phaseaddp'+str(self.fitLew)] = pf[0] %(2*np.pi)
		if (self.fitLew) : self.Lewp = pf[1]
		#--------------------------------------------------
		if (self.fitLew) : np.save(self.outdir+'Lewp_'+str(self.nsigma)+'_m'+str(self.fitLew)+'.npy', self.Lewp)
		np.save(self.outdir+'Phaseaddp'+str(self.fitLew)+'_'+str(self.nsigma)+'.npy', self.__dict__['Phaseaddp'+str(self.fitLew)])
		if (self.verbose) : 
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
		outdir = jp.Mkdir(self.outdir + 'Plot/')
		#--------------------------------------------------
		nbl = np.arange(len(self.antarray[0].visorder))
		bl = self.antarray[0].Blorder.blorder[self.antarray[0].visorder]
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
		freq = self.antarray[0].Ant.freq
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
		bcast = (freq, func, color, outdir, 'amp', self.nsigma)
		if (Nprocess <= 1) : 
			iterable = [(0,len(strbl)), send, bcast]
			_DoMultiprocess_Plot( iterable )
		else : 
			pool = jp.PoolFor(0, len(strbl), Nprocess)
			pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotDeff( self, outdir, strbl, Nprocess, dyDeff ) : 
		freq = self.antarray[0].Ant.freq
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
		antform = jp.npfmt(self.antarray[0].Ant.antform)[0]
		bcast = (freq, func, color, outdir, 'deff', self.nsigma, dyDeff, self.fwhm2sigmafactor, antform)
		if (Nprocess <= 1) : 
			iterable = [(0,len(strbl)), send, bcast]
			_DoMultiprocess_Plot( iterable )
		else : 
			pool = jp.PoolFor(0, len(strbl), Nprocess)
			pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotLew( self, outdir, strbl, Nprocess, dyLew ) : 
		freq = self.antarray[0].Ant.freq
		func, Lew, color = ['FitPhase', 'FitVis'], [], ['b','r']
		try    : Lew.append(self.Lewp.T)  # (Nv,Nf)
		except : Lew.append(None)
		try    : Lew.append(self.Lewv.T)
		except : Lew.append(None)
		if (Lew[0] is None and Lew[1] is None) : 
			if (self.verbose) : print '    NOT exist  self.Lewp, self.Lewv'
			return
		if (self.verbose) : print '    Plotting    Lew     ......'
		bl =self.antarray[0].Blorder.baseline[self.antarray[0].visorder]
		send = tuple(Lew) + (strbl,)
		bcast = (freq, func, color, outdir, 'lew', self.nsigma, self.fitLew, bl, dyLew)
		if (Nprocess <= 1) : 
			iterable = [(0,len(strbl)), send, bcast]
			_DoMultiprocess_Plot( iterable )
		else : 
			pool = jp.PoolFor(0, len(strbl), Nprocess)
			pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotPhaseadd( self, outdir, strbl, Nprocess ) : 
		freq = self.antarray[0].Ant.freq
		func, phaseadd, color = ['FitPhase-m1','FitPhase-m2','FitVis','FitExp'], [], ['b','r','c','k']
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
		bl =self.antarray[0].Blorder.baseline[self.antarray[0].visorder]
		send = tuple(phaseadd) + (strbl,)
		bcast = (freq, func, color, outdir, 'phaseadd', self.nsigma, self.fitLew, bl)
		if (Nprocess <= 1) : 
			iterable = [(0,len(strbl)), send, bcast]
			_DoMultiprocess_Plot( iterable )
		else : 
			pool = jp.PoolFor(0, len(strbl), Nprocess)
			pool.map_async(_DoMultiprocess_Plot, send, bcast)



	def _PlotFringe( self, array, label, color, lw=None ) : 
		'''
		array: 
			vis to be plotted
			(1) ndarray
			(2) [ndarray, ndarray, ...], in this case, label, color, lw should has the same number

		label: 
			plt.plot(x, y, label=label)
			label must be str pair [label1, label2], 1=>real, 2=>imag

		color:
			color must be str pair [color1, color2], 1=>real, 2=>imag

		lw:
			linewidth, must be int, one int for real+imag
			can be lw=None
		'''
		# Plot, longest East-West baseline, center frequency
		if (self.timem[-1]-self.timem[0] <= 60) : xmajor = 5
		else : xmajor = 10
		#--------------------------------------------------
		istype = jp.IsType()
		if (not (istype.islist(array) or istype.istuple(array))) : 
			if (lw is None) : lw = 1 
			array, label, color, lw= [array], [label], [color], [lw]
		else : 
			if (lw is None) : lw = [1,1] 
		#--------------------------------------------------
		if (self.plotnv is None) : 
			bl = self.antarray[0].Blorder.baseline[self.antarray[0].visorder]
			bl = abs(bl[:,0]) # longest East-West baseline
			nv = np.where(bl==bl.max())[0][0]
		else : nv = self.plotnv
		freq=self.antarray[0].Ant.freq[self.antarray[0].freqorder]
		if (self.plotfreq is None) : plotfreq = freq.mean()
		else : plotfreq = self.plotfreq
		nvbl = np.arange(len(self.antarray[0].Blorder.blorder))[self.antarray[0].visorder][nv]
		bl = self.antarray[0].Blorder.Order2Bl(nvbl)
		strbll = '=(%.3f, %.3f)' % tuple(self.antarray[0].Blorder.baseline[nvbl][:2])
		strbl = str(bl[0])+'-'+str(bl[1])
		nf = abs(freq-plotfreq)
		nf = np.where(nf==nf.min())[0][0]
		freqstr = str(int(round(freq[nf])))
		#--------------------------------------------------
		# Plot fringe
		vmax, vmin = [], []
		if (self.vistype == 'cross') : plt.figure(figsize=(17,6))
		for i in xrange(len(array)) : 
			vmax.append( abs(array[i][:,nf,nv]).max()*1.05 )
			vmin.append( abs(array[i][:,nf,nv]).min()*0.95 )
			scale = jp.SciNot(vmax[-1])[1]-2
			vmax[-1] /= 10.**scale
			vmin[-1] /= 10.**scale
			if (self.vistype == 'cross') : 
				plt.subplot(1,2,1)
				plt.plot(self.timem, array[i].real[:,nf,nv]/10.**scale, color=color[i][0], lw=lw[i], label=label[i][0])
			else : 
				plt.plot(self.timem, array[i][:,nf,nv]/10.**scale, color=color[i][0], lw=lw[i], label=label[i][0])
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
			plt.subplot(1,2,2)
			plt.plot(self.timem, array[i].imag[:,nf,nv]/10.**scale, color=color[i][1], lw=lw[i], label=label[i][1])
			if (i == len(array)-1) : 
				plt.legend()
				plt.xlabel('time [min]', size=16)
				plt.xlim(self.timem.min(), self.timem.max())
				plt_axes('x', 'both', [xmajor,1])
				plt.ylabel('[A.U.]', size=16)
				plt.ylim(-vmax, vmax)
		if (self.vistype == 'cross') : plt.suptitle('Fringe of baseline='+strbl+strbll+' @ '+freqstr+'MHz', fontsize=20)
		else : 
			strbl = strbl[:strbl.find('-')]
			plt.suptitle('Auto-correlation of channel='+strbl+' @ '+freqstr+'MHz', fontsize=16)
		plt.savefig(self.outdir+'vis_'+strbl+'_'+freqstr+'MHz_'+str(self.nsigma)+'.png')
		plt.close()
		if (self.vistype != 'cross') : return
		#--------------------------------------------------
		# Plot phase
		label = ['Phase', 'Smoothed']
		markersize = [4, 3]
		colorls = ['bo', 'ro']
		for i in xrange(len(array)) : 
			phase = np.angle(array[i][:,nf,nv])
			phase = (phase*180/np.pi) %360
			plt.plot(self.timem, phase, colorls[i], markersize=markersize[i], label=label[i])
		plt.legend(fontsize=10)
		plt.xlabel('time [min]', size=16)
		plt.xlim(self.timem.min(), self.timem.max())
		plt_axes('x', 'both', [xmajor,1])
		plt.ylim(0, 360)
		plt.ylabel('[degree]', size=16)
		plt.title('Total phase of baseline='+strbl+strbll+' @ '+freqstr+'MHz', size=13)
		plt.savefig(self.outdir+'phase_'+strbl+'_'+freqstr+'MHz_'+str(self.nsigma)+'.png')
		plt.close()



