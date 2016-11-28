import jizhipy as jp
from jizhipy.Plot import *
from AntArray import *
from Masking import *
##################################################



class CaliTemp( object ) : 


	def __init__( self, caliphasecross=None, caliphaseauto=None, verbose=True, outdir=None ) : 
		'''
		caliphasecross:
			caliphasecross.vis: with masking, caligain
		'''
		self.starttime = jp.Time(1)
		self.verbose = verbose
		if (self.verbose) : print '\n'
		if (caliphasecross is not None) : self.caliphasecross = caliphasecross
		if (caliphaseauto is not None) : self.caliphaseauto = caliphaseauto
		self.outdir = jp.Outdir((None,'file'), (0,'file'))
		if (outdir is not None) : self.outdir = jp.Mkdir(self.outdir+outdir)



	def Calibrator( self, fluxdensity=None, sourcename=None ) : 
		'''
		Just need fluxdensity or sourcename (one of two)
		Use fluxdensity first, if not, then use sourcename

		fluxdensity:
			in Jy, the same shape as AntArray.Ant.freq

		sourcename:
			Format of the sourcename, CygA as an example:
				sourcename = 'CygA' or 'Cygnus_A'
		'''
		if (fluxdensity is not None) : self.flux = fluxdensity
		elif (sourcename is not None) : 
			bs = jp.BrightSource()
			self.flux = bs.FluxDensity(sourcename, self.caliphasecross.antarray[0].Ant.freq)



	def Temp( self, kDeff=1 ) : 
		'''
		kDeff:
		We use caliphasecross.FitBeam to fit Deff. If think the result not good, use kDeff*Deff to instead

		return:
		self.calitempcross:
			Basing on caliphasecross: whether with caligaint
		'''
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'CaliTemp.Temp: start @', starttime
		# Temperature of the source in K
		Deff = (300/self.caliphasecross.antarray[0].Ant.freq[:,None] /self.caliphasecross.Sigmab /self.caliphasecross.fwhm2sigmafactor).mean(1)
		n1, n2 = Deff.size/5, Deff.size/5*4
		x = np.arange(n2-n1)
		def func(x, p) : return p[0]
		Deff = jp.Leastsq(func, x, Deff[n1:n2], [0,1])[0]
		Ts = jp.Jy2K(self.flux, Deff*kDeff)[:,None]
		#--------------------------------------------------
		vistype = self.caliphasecross.antarray[0].vistype[:-1]
		self.calitempcross = Ts /self.caliphasecross.Ampb
		print 'Ts.shape =', Ts.shape
		print 'self.caliphasecross.Ampb.shape =', self.caliphasecross.Ampb.shape
		print 'self.calitempcross.shape', self.calitempcross.shape
		np.save(self.outdir+'calitempcross', self.calitempcross)
		print 'self.caliphasecross.caligain[0].norm.shape', self.caliphasecross.caligain[0].norm.shape
		#--------------------------------------------------
		self.Tnoisesigmacross = self.caliphasecross.caligain[0].norm *self.calitempcross
		print 'self.Tnoisesigmacross.shape', self.Tnoisesigmacross.shape
		print 'self.caliphaseauto.caligain[0].norm.shape', self.caliphaseauto.caligain[0].norm.shape
		self.calitempauto = self.caliphaseauto.caligain[0].norm * 2**0.5 *self.calitempcross 
		print 'calitempauto.shape', calitempauto.shape
		np.save(self.outdir+'calitempauto', self.calitempauto)
		np.save(self.outdir+'Tnoisesigmacross', self.Tnoisesigmacross)
		print '    calitempcross.npy, calitempauto.npy, Tnoisesigmacross.npy  saved'
		#--------------------------------------------------
		if (self.verbose) : 
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliTemp.Temp:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'



	def Plot( self, plotfreq=None, plotnv=None, yaxes=None, legendsize=None, legendloc=None ) : 
		'''
		plotfreq:
			in MHz
			==None: middle frequency

		plotnv:
			Used to Cross-correlation, NOT Auto
			caliphasecross.antarray[0].visorder[plotnv]
			plotnv is the number/index/order of visorder, NOT the number/index/order of AntArray.Blorder.blorder
			==None: longest East-West baseline

		legendsize, legendloc:
			Can be one value: legendsize=10, for both auto and cross
			OR pair: legendsize=[10,12], 10 for auto, 12 for cross
		'''
		if (plotnv is None) : # longest East-West baseline
			bl = self.caliphasecross.antarray[0].Blorder.baseline[self.caliphasecross.antarray[0].visorder]
			bl = abs(bl[:,0]) # longest East-West baseline
			nvc = np.where(bl==bl.max())[0][0]  #@#@
		else : nvc = plotnv  #@#@
		nvc = jp.npfmt(nvc)
		bli = self.caliphasecross.antarray[0].Blorder.baseline[self.caliphasecross.antarray[0].visorder][nvc]
		blo = self.caliphasecross.antarray[0].Blorder.blorder[self.caliphasecross.antarray[0].visorder][nvc]
		#--------------------------------------------------
		chan, nva = [], []
		bla = self.caliphaseauto.antarray[0].Blorder.blorder[self.caliphaseauto.antarray[0].visorder][:,0]
		for i in xrange(nvc.size) : 
			chan1, chan2 = blo[i]
			if (chan1 not in chan) : 
				chan.append(chan1)
				nva.append(np.where(bla==chan1)[0][0])
			if (chan2 not in chan) : 
				chan.append(chan2)
				nva2 = np.where(bla==chan2)[0][0]
		nva = np.array(nva)
		#--------------------------------------------------
		freq = self.caliphasecross.antarray[0].Ant.freq
		if (plotfreq is None) : plotfreq = freq.mean()
		freqstr = str(int(round(plotfreq)))  #@#@
		nf = abs(freq-plotfreq)
		nf = np.where(nf==nf.min())[0][0]  #@#@
		#--------------------------------------------------
		viscross = self.caliphasecross.vis[:,nf,nvc] *self.calitempcross[nf,nvc]
		print viscross.shape
		print self.caliphaseauto.vis.shape
		print nf
		print nva, nva.shape
		visauto = self.caliphaseauto.vis[:,nf,nva]*self.calitempauto[nf,nva]
		print visauto.shape
		timem = self.caliphasecross.timem
		legendsize, legendloc = jp.npfmt(legendsize), jp.npfmt(legendloc)
		if (legendsize.size == 1) : 
			legendsize = np.append(legendsize, legendsize)
			legendloc = np.append(legendloc, legendloc)
		plt.figure(figsize=(17,6))
		#--------------------------------------------------
		# Plot all auto
		plt.subplot(1,2,1)
		color = plt_color(visauto.shape[-1])
		ymin, ymax = visauto.min(), visauto.max()
		dy = ymax-ymin
		ymin = ymin-0.1*dy if(ymin>0.1*dy)else 0
		ymax += 0.1*dy
		for i in xrange(visauto.shape[-1]) : 
			plt.plot(timem, visauto[:,i], color=color[i], label=str(chan[i]))
		plt.legend(loc=legendloc[1], fontsize=legendsize[1])
		plt.xlim(timem.min(), timem.max())
		plt.ylim(ymin, ymax)
		yaxes
		plt.xlabel(r't [min]', size=16)
		plt.ylabel(r'T [K]', size=16)
		plt.title('Auto-correlations', size=16)
		#--------------------------------------------------
		plt.subplot(1,2,2)
		color = plt_color(2*viscross.shape[-1])
		ymax = abs(viscross).max()
		for i in xrange(viscross.shape[-1]) : 
			strbli = '(%.3f,%.3f)' % tuple(bli[:2,i])
			strblo = '%i-%i' % tuple(blo[:,i])
			plt.plot(timem, viscross[:,i].real, color=color[i], label=strblo+', real')
			plt.plot(timem, viscross[:,i].imag, color=color[-i-1], label=strblo+', imag')
		plt.legend(loc=legendloc[1], fontsize=legendsize[1])
		plt.xlim(timem.min(), timem.max())
		plt.ylim(-1.05*ymax, 1.05*ymin)
		yaxes
		plt.xlabel(r't [min]', size=16)
		plt.ylabel(r'T [K]', size=16)
		if (nvc.size == 1) : 
			plt.title('Cross-correlation of baseline='+strblo+'='+strbli, size=16)
		else : plt.title('Cross-correlations', size=16)
		plt.suptitle(r'Calibrating the amplitude to Kelvin @ '+str(self.freq)+'MHz', size=20)
		if (nvc.size == 1) : 
			plt.savefig(self.outdir+'vis_'+strblo+'_'+str(self.freq)+'MHz.png')
		else : plt.savefig(self.outdir+'vis_'+str(self.freq)+'MHz.png')
		plt.close()



