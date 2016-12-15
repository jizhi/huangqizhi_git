import jizhipy as jp
from jizhipy.Plot import *
from AntArray import *
from Masking import *
##################################################



class CaliTemp( object ) : 



	def __init__( self, caliphasecross, caliphaseauto=None, outdir='' ) : 
		'''
		Must give caliphasecross.
		If also want to calibrate auto, give caliphaseauto
		'''
		self.verbose = caliphasecross.verbose
		if (self.verbose) : print '-------------------- CaliTemp --------------------\n'
		self.starttime = jp.Time(1)
		self.caliphasecross = caliphasecross
		if (caliphaseauto is not None) : self.caliphaseauto = caliphaseauto
		self.outdir = jp.Outdir((None,'file'), (0,'file'))
		if (outdir is not None) : self.outdir = jp.Mkdir(self.outdir+outdir)
		try : self.outdir = jp.Mkdir(self.outdir + outdir)
		except : pass





	def CalibratorFlux( self, fluxdensity=None, sourcename=None ):
		'''
		Just need fluxdensity or sourcename (one of two)
		Use fluxdensity first, if not, then use sourcename

		fluxdensity:
			in Jy, the same shape as Antarray.freqorder

		sourcename:
			Format of the sourcename, CygA as an example:
				sourcename = 'CygA' or 'Cygnus_A'
			Then use jizhipy.BrightSource() to get fluxdensity
		'''
		if (fluxdensity is not None) : self.flux = fluxdensity
		elif (sourcename is not None) : 
			freq = self.caliphasecross.antarray.Ant.freq[self.caliphasecross.antarray.freqorder]
			bs = jp.BrightSource()
			self.flux = bs.FluxDensity(sourcename, freq)





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
			print 'CaliTemp.Temp:  start @', starttime
		# Temperature of the source in K
		freq = self.caliphasecross.antarray.Ant.freq[self.caliphasecross.antarray.freqorder]
		try : 
			Deff = (300/freq[:,None] /self.caliphasecross.Sigmab /self.caliphasecross.fwhm2sigmafactor).mean(1)  # Assume all cross have the same Deff
			n1, n2 = Deff.size/5, Deff.size/5*4
			x = np.arange(n2-n1)
			def func(x, p) : return p[0]
			Deff = jp.Leastsq(func, x, Deff[n1:n2], [0])[0]
		except: Deff= self.caliphasecross.antarray.Ant.antform*0.9
		Ts = jp.Jy2K(self.flux, freq, Deff*kDeff)[:,None]
		#--------------------------------------------------
		vistype = self.caliphasecross.antarray.vistype[:-1]
		self.Tcross = Ts /self.caliphasecross.Ampb
		np.save(self.outdir+'Tcross.npy', self.Tcross)
		#--------------------------------------------------
		self.Tnoisesigmacross = self.caliphasecross.caligain.norm *self.Tcross
		np.save(self.outdir+'Tnoisesigmacross', self.Tnoisesigmacross)
		try : 
			self.Tauto = self.caliphaseauto.caligain.norm * 2**0.5 *self.Tcross 
			np.save(self.outdir+'Tauto', self.Tauto)
			printstr = '    calitempcross.npy, calitempauto.npy, Tnoisesigmacross.npy  saved to '+self.outdir
		except : printstr = '    calitempcross.npy, Tnoisesigmacross.npy  saved to '+self.outdir
		#--------------------------------------------------
		if (self.verbose) : 
			print printstr
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'CaliTemp.Temp:   end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'





	def Plot( self, nf=None, nv=None, yaxes=None, legendsize=None, legendloc=None, blankline=False ) : 
		'''
		nf, nv:
			Both single int
			subplot(1,2,1) plots 1 cross
			subplot(1,2,2) plots 2 autos

		yaxes:
			[yaxescross, yaxesauto]
			If not plot auto, then yaxes=yaxescross, or yaxes=[yaxescross]
			Any element can be None

		legendsize, legendloc:
			Can be one: then copy to two
 			Can be two: [legendsizecross, legendsizeauto], .....
			Any element can be None

		blankline:
			==True: print at the end
		'''
		nf, strfreq, nv, strbloc, strblic, blnamec = self.caliphasecross.antarray.Plotnfnv(nf, nv, False)
		bloc = np.sort(jp.npfmt(strbloc.split('-'), int))
		try : 
			strbloa, strblia, blnamea = [], [], []
			bloa = self.caliphaseauto.antarray.Blorder.blorder[self.caliphaseauto.antarray.visorder][:,0]
			na = np.concatenate(np.where(bloa==bloc[0]) + np.where(bloa==bloc[1]))
			for i in xrange(len(na)) : 
				nf1, strfreq1, na1, strbloa1, strblia1, blnamea1 = self.caliphase.antarray.Plotnfnv(nf, na[i], False)
				strbloa.append(strbloa1)
				strblia.append(strblia1)
				blnamea.append(blnamea1)
		except : pass
		#--------------------------------------------------
		istype = jp.IsType()
		if (yaxes is None) : yaxes = [None, None]
		elif (not istype.islist(yaxes) and not istype.istuple(yaxes)) : yaxes = [yaxes, None]
		#--------------------------------------------------
		if (legendsize is None) : legendsize = [None, None]
		elif (not istype.islist(legendsize) and not istype.istuple(legendsize)) : legendsize = [legendsize, legendsize]
		#--------------------------------------------------
		if (legendloc is None) : legendloc = [None, None]
		elif (not istype.islist(legendloc) and not istype.istuple(legendloc)) : legendloc = [legendloc, legendloc]
		#--------------------------------------------------
		viscross = self.caliphasecross.antarray.vis[:,nf,nv] * self.Tcross[nf,nv]
		maskcross = self.caliphasecross.masking.mask[:,nf,nv]
		viscross = np.ma.MaskedArray(viscross, maskcross)
		try : 
			visauto =self.caliphaseauto.antarray.vis[:,nf,na] * self.Tauto[nf,nv]
			try : 
				maskauto = self.caliphaseauto.masking.mask[:,nf,na]
				viscauto = np.ma.MaskedArray(viscauto, maskcross)
			except : pass
		except : pass
		#--------------------------------------------------
		plt.figure(figsize=(17,6))
		plt.subplot(1,2,1)
		timem = self.caliphasecross.antarray.timem
		plt.plot(timem, viscross.real, 'b-', label='data.real')
		plt.plot(timem, viscross.imag, 'r-', label='data.imag')
		try : plt.legend(loc=legendloc[0], fontsize=legendsize[0])
		except : pass
		plt.xlabel('time [min]', size=16)
		plt.xlim(timem[0], timem[-1])
		vmax = abs(viscross).max()
		plt.ylim(-vmax, vmax)
		try : yaxes[0]
		except : pass
		plt.ylabel(r'$T_b$ [K]', size=16)
		plt.title('Cross-correlation', size=16)
		#--------------------------------------------------
		try : 
			plt.subplot(1,2,2)
			timem = self.caliphaseauto.antarray.timem
			color = ['b-', 'r-']
			for i in xrange(len(na)) : 
				plt.plot(timem, visauto[:,i], color[i], label='chan '+strbloa[i])
			try: plt.legend(loc=legendloc[1],fontsize=legendsize[1])
			except : pass
			plt.xlabel('time [min]', size=16)
			plt.xlim(timem[0], timem[-1])
			try : yaxes[1]
			except : pass
			plt.ylabel(r'$T_b$ [K]', size=16)
			plt.title('Auto-correlation', size=16)
		except : pass
		plt.suptitle('CaliTemp, '+strfreq+', '+blnamec+strbloc+strblic, size=16)
		plt.savefig(self.outdir+'CaliTemp_'+strfreq+'_'+strbloc+'.png')
		plt.close()
		if (self.verbose) : print '    Plotting CaliTemp_'+strfreq+'_'+strbloc+'.png'
		if (blankline) : print



