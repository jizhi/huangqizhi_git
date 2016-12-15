import jizhipy as jp
from jizhipy.Plot import *
import healpy as hp
import pyfits
import numpy as np
import os
from AntArray import *





def _DoMultiprocess_Vis( iterable ) : 
	Nmax = iterable[0][1]
	timelist = iterable[1]
	tmap2vis = Tmap2Vis(verbose=False)
	RAstart, Bmap, freq, tmap2vis.hpixCeleXYZ, tmap2vis.blXYZ, tmap2vis.Tmap, tmap2vis.verbose = iterable[2]
	coordtrans = jp.CoordTrans()
	if (tmap2vis.verbose) : 
		progressbar = jp.ProgressBar('        Completed', Nmax, False)
		for i in xrange(Nmax-len(timelist)): progressbar.Progress()
	#--------------------------------------------------
	vis = []
	for j in xrange(timelist.size) : 
		if (tmap2vis.verbose) : progressbar.Progress()
		dRA = timelist[j]/3600.*15
		RAnow = RAstart + dRA
		BmapComplex = tmap2vis.BeamMapComplex(RAnow, Bmap, freq)  # (1, Nv, 12*nside**2)
		visi = BmapComplex * tmap2vis.Tmap[None,None,:]
		BmapComplex = 0 #@
		visi = visi.sum(-1)  # Sum all directions
		vis.append( visi[None,:] )
		#--------------------------------------------------
	vis = np.concatenate(vis, 0)  # (Nt, 1, Nv)
	return vis





class Tmap2Vis( object ) : 



	def __init__( self, nside=None, Nprocess=None, verbose=True, outdir='' ) :
		'''
		nside should be larger than Tmap.nside
		nside >= 4*Tmap.nside
		'''
		self.starttime = jp.Time(1)
		self.verbose, self.Nprocess = bool(verbose), Nprocess
		if (self.verbose) : print '-------------------- Tmap2Vis --------------------\n'
		self.Nprocess = jp.NprocessCPU(Nprocess, verbose)[0]
		self.outdir = jp.Outdir((None,'file'), (0,'file'))
		self.nside = nside
		try : self.outdir = jp.Mkdir(self.outdir + outdir)
		except : pass



################### HealpixMap START #####################



	def HealpixMap( self, inpath, which=(1,0), ordering=None, coordsys=None, unit=None, nside=None, dtype=None ) : 
		'''
		Input map must be healpix map

		healmap:
			(1) Path of healpix map with ".fits": 
				hpmap = hp.read_map(inpath, hdu=which[0], field=which[1], nest=nest, dtype=dtype, verbose=verbose)
				which = (hdu, field)

			(2) Path of healpix map with ".hdf5": 
				hpmap = h5py.open(inpath)[which]
				which is the key (str)

			(3) np.ndarray with healpix pixelization (size=12*2**n)

		self.Tmap
		'''
		if (self.verbose) : print 'Tmap2Vis.HealpixMap'
		self.ordering, self.nest = self._OrderingCheck(None, ordering)[2:]
		self.coordsys = self._CoordsysCheck(None, coordsys)[1]
		if (nside is not None) : 
			self.nside = self._NsideCheck(nside)
		self.unit = self._UnitCheck(None, unit)[1]
		if (dtype is not None) : 
			try    : dtype = np.dtype(dtype)
			except : dtype = None
		self.dtype = dtype
		#--------------------------------------------------

		istype = jp.IsType()
		if (not istype.isstr(inpath)) : 
			self.Tmap = np.array(inpath)
			self.nsidein = hp.get_nside(self.Tmap)
			self.dtypein = self.Tmap.dtype
			self.Tmap = np.array(self.Tmap, self.dtype)
			self.dtype = self.Tmap.dtype
			self.inpath = None
			self.orderingin = None
			self.coordsysin = None
			self.unitin = None
			self.freqin = None
			printstrNone = ['orderingin', 'coordsysin', 'unitin', 'freqin']
		#--------------------------------------------------

		elif ('.fits' == str(inpath).lower()[-5:]) : 
			printstrNone = []
			self.inpath =os.path.abspath(os.path.expanduser(inpath))
			fo = pyfits.open(self.inpath)
			hdr = fo[0].header
			for i in xrange(1, which[0]+1): hdr.update(fo[i].header)
			fo = fo[which[0]]
			self.dtypein = np.dtype(fo.data.dtype[which[1]].name)
			self.Tmap = np.array(fo.data[fo.data.names[which[1]]], self.dtype)  # Read data directly
			self.dtype = self.Tmap.dtype
			self.nsidein = hp.get_nside(self.Tmap)
			#--------------------------------------------------
			try : 
				orderingin = str(hdr['ORDERING']).upper()
				orderingin, nestin = self._OrderingCheck(orderingin)[:2]
			except : 
				printstrNone.append('orderingin')
				orderingin, nestin = None, None
			self.orderingin, self.nestin = orderingin, nestin
			#--------------------------------------------------
			try : 
				try : coordsysin = str(hdr['COORDSYS']).lower()
				except : coordsysin = str(hdr['SKYCOORD']).lower()
				coordsysin = self._CoordsysCheck(coordsysin)[0]
			except : 
				printstrNone.append('coordsysin')
				coordsysin = None
			self.coordsysin = coordsysin
			#--------------------------------------------------
			unitin = 'None'
			try : 
				unitin = str(hdr['TUNIT'+str(which[1]+1)])
				unitin = self._UnitCheck(unitin)[0]
			except : 
				print "    Warning: unitin = "+unitin+" not in ['K', 'mK', 'uK']"
				printstrNone.append('unitin')
				unitin = None
			self.unitin = unitin
			#--------------------------------------------------
			printstr = ''
			try : 
				n = hdr.keys().index('FREQ')
				freqin, com = str(hdr['FREQ']).lower(), hdr.comments[n].lower()
				if   ('ghz' in freqin or 'ghz' in com) : freqin = 1e3 * float(freqin.split('ghz')[0])
				elif ('mhz' in freqin or 'mhz' in com) : freqin = float(freqin.split('mhz')[0])
				elif ( 'hz' in freqin or  'hz' in com) : freqin = float(freqin.split('hz')[0])
				else : printstr = 'FREQ : '+freqin+' / '+com+", NOT with unit ['Hz', 'MHz', 'GHz']"
			except : 
				printstrNone.append('freqin')
				freqin = None
			self.freqin = freqin
		#--------------------------------------------------

		if (len(printstrNone) > 0) : 
			printstr = ''
			for ps in printstrNone : printstr += ps + ', '
		if (printstr != '') : print '    Warning: '+printstr[:-2]+' = None'
		if (self.orderingin is None) : "    Warning: orderingin = None, can NOT change to ordering ('NESTED'<=>'RING'), nside (ud_grade), coordsys ('galactic'<=>'equatorial') in __init__()"
		#--------------------------------------------------

		if (self.orderingin is not None) : 
			# coordsys
			if (self.coordsysin is None and self.coordsys is not None) : 
				print "    Warning: coordsysin = None, coordsys = '"+self.coordsys+"', can NOT change"
				self.coordsys = None
			elif (self.coordsysin is None and self.coordsys is None) : pass
			elif (self.coordsysin is not None and self.coordsys is None) : self.coordsys = self.coordsysin
			elif (self.coordsysin != self.coordsys) : self.CoordsysConvert() 
			#--------------------------------------------------

			# ordering
			if (self.ordering is None) : self.ordering, self.nest = self.orderingin, self.nestin
			if (self.orderingin != self.ordering) : self.OrderingConvert()
			#--------------------------------------------------

			# nside
			if (self.nside is None) : self.nside = self.nsidein
			if (self.nsidein != self.nside): self.NsideConvert(None, self.ordering)
		#--------------------------------------------------

		if (self.unitin is None) : self.unit = None
		elif (self.unit is None) : self.unit = self.unitin
		elif (self.unitin != self.unit) : self.UnitConvert()
		if (self.verbose) : print
		# Tmap2Vis.HealpixMap END





	def _CoordsysCheck( self, coordsysin=None, coordsysout=None): 
		''' '''
		doraise = []
		if (coordsysin is None) : 
			try : coordsysin = self.coordsysin
			except : pass
		else : 
			coordsysin = str(coordsysin)
			if (coordsysin.lower() not in ['galactic','equatorial']) : doraise.append('in')
			else : coordsysin = coordsysin.lower()
		#--------------------------------------------------
		if (coordsysout is None) : 
			try : coordsysout = self.coordsys
			except : pass
		else : 
			coordsysout = str(coordsysout)
			if (coordsysout.lower() not in ['galactic','equatorial']) : doraise.append('out')
			else : coordsysout = coordsysout.lower()
		#--------------------------------------------------
		printstr = ''
		if ('in' in doraise) : printstr += "coordsysin = '"+coordsysin+"', "
		if ('out' in doraise) : printstr += "coordsysout = '"+coordsysout+"', "
		if (printstr != '') : jp.Raise(Exception, printstr[:-2]+" not in ['galactic', 'equatorial']")
		else : return [coordsysin, coordsysout]





	def _OrderingCheck( self, orderingin=None, orderingout=None):
		''' '''
		doraise = []
		if (orderingin is None) : 
			try : orderingin = self.orderingin
			except : pass
		else : 
			orderingin = str(orderingin)
			if (orderingin.upper() not in ['RING', 'NESTED', 'NEST']) : doraise.append('in')
			else : 
				orderingin = orderingin.upper()
				if ('NEST' in orderingin) : orderingin = 'NESTED'
		#--------------------------------------------------
		if (orderingout is None) : 
			try : orderingout = self.ordering
			except : pass
		else : 
			orderingout = str(orderingout)
			if (orderingout.upper() not in ['RING', 'NESTED', 'NEST']) : doraise.append('out')
			else : 
				orderingout = orderingout.upper()
				if ('NEST' in orderingout) : orderingout = 'NESTED'
		#--------------------------------------------------
		printstr = ''
		if ('in' in doraise) : printstr += "orderingin = '"+orderingin+"', "
		if ('out' in doraise) : printstr += "orderingout = '"+orderingout+"', "
		if (printstr != '') : jp.Raise(Exception, printstr[:-2]+" not in ['RING', 'NESTED']")
		#--------------------------------------------------
		if (orderingin is None) : nestin = None
		elif (orderingin == 'RING') : nestin = False
		else : nestin = True
		if (orderingout is None) : nestout = None
		elif (orderingout == 'RING') : nestout = False
		else : nestout = True
		return [orderingin, nestin, orderingout, nestout]





	def _NsideCheck( self, nsideout=None ) : 
		if (nsideout is None) : 
			try : nsideout = self.nside
			except : pass
		else : 
			n = np.log(nsideout) / np.log(2)
			reset = True if(n != int(n))else False
			if (reset) : print '   Warning: nsideout = 2**'+('%.3f' % n)+' != 2**int, ',
			n = int(round(n))
			nsideout = 2**n
			if (reset) : print 'reset nsideout = 2**'+str(n)+' = '+str(nsideout)
		return nsideout





	def _UnitCheck( self, unitin=None, unitout=None ) : 
		''' return [unitin, unitout, k_in2out] '''
		doraise = []
		if (unitin is None) : 
			try : unitin = self.unitin
			except : pass
		else : 
			unitin = str(unitin)
			if   (unitin.lower() ==  'k') : unitin =  'K'
			elif (unitin.lower() == 'mk') : unitin = 'mK'
			elif (unitin.lower() == 'uk') : unitin = 'uK'
			else : doraise.append('in')
		#--------------------------------------------------
		if (unitout is None) : 
			try : unitout = self.unit
			except : pass
		else : 
			unitout = str(unitout)
			if   (unitout.lower() ==  'k') : unitout =  'K'
			elif (unitout.lower() == 'mk') : unitout = 'mK'
			elif (unitout.lower() == 'uk') : unitout = 'uK'
			else : doraise.append('out')
		#--------------------------------------------------
		printstr = ''
		if ('in' in doraise): printstr+="unitin = '"+unitinin+"', "
		if('out' in doraise): printstr+="unitout = '"+unitout+"', "
		if (printstr != '') : jp.Raise(Exception, printstr[:-2]+" not in ['K', 'mK', 'uK']")
		#--------------------------------------------------
		if (unitin is None or unitout is None): k_in2out = 1
		elif (unitin== 'K' and unitout== 'K') : k_in2out = 1
		elif (unitin== 'K' and unitout=='mK') : k_in2out = 1e3
		elif (unitin== 'K' and unitout=='uK') : k_in2out = 1e6
		elif (unitin=='mK' and unitout== 'K') : k_in2out = 1e-3
		elif (unitin=='mK' and unitout=='mK') : k_in2out = 1
		elif (unitin=='mK' and unitout=='uK') : k_in2out = 1e3
		elif (unitin=='uK' and unitout== 'K') : k_in2out = 1e-6
		elif (unitin=='uK' and unitout=='mK') : k_in2out = 1e-3
		elif (unitin=='uK' and unitout=='uK') : k_in2out = 1
		return [unitin, unitout, k_in2out]





	def CoordsysConvert( self, coordsysin=None, coordsysout=None, orderingin=None, inmap=None ) : 
		''' '''
		coordsysin, coordsysout = self._CoordsysCheck(coordsysin, coordsysout)
		orderingin = self._OrderingCheck(orderingin)[0]
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : strfont = ''
		else : strfont = '    '
		if (self.verbose) : 
			print strfont+'Tmap2Vis.CoordsysConvert'
			print strfont+strfont+"coordsysin='"+coordsysin+"', coordsysout='"+coordsysout+"', ordering='"+orderingin+"'"
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : print
		coordtrans = jp.CoordTrans(self.Nprocess)
		if (inmap is None) : 
			self.Tmap = coordtrans.CelestialHeal(self.Tmap, coordsysin, coordsysout, orderingin)
		else : return coordtrans.CelestialHeal(inmap, coordsysin, coordsysout, orderingin)




	def OrderingConvert( self, orderingin=None, orderingout=None, inmap=None ) : 
		''' '''
		orderingin, nestin, orderingout = self._OrderingCheck(orderingin, orderingout)[:3]
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : strfont = ''
		else : strfont = '    '
		if (self.verbose) : 
			print strfont+'Tmap2Vis.OrderingConvert'
			print strfont+strfont+"orderingin='"+orderingin+"', orderingout='"+orderingout+"'"
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : print
		if (inmap is None) : self.Tmap = hp.ud_grade(self.Tmap, self.nsidein, order_in=orderingin, order_out=orderingout)
		else : return hp.ud_grade(inmap, hp.get_nside(inmap), order_in=orderingin, order_out=orderingout)





	def NsideConvert( self, nsideout=None, orderingin=None, inmap=None ) :
		''' '''
		nsideout = self._NsideCheck(nsideout)
		orderingin = self._OrderingCheck(orderingin)[0]
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : strfont = ''
		else : strfont = '    '
		if (self.verbose) : 
			print strfont+'Tmap2Vis.NsideConvert'
		nsideinstr = ''
		if (inmap is None) : 
			print strfont+strfont+"nsidein="+str(hp.get_nside(self.Tmap))+", nsideout="+str(nsideout)+", orderingin='"+orderingin+"'"
			if(self.verbose and jp.SysFrame(0,2)[-1][-2]==''): print
			self.Tmap = hp.ud_grade(self.Tmap, nsideout, order_in=orderingin)
		else : 
			print strfont+strfont+"nsidein="+str(hp.get_nside(inmap))+", nsideout="+str(nsideout)+", orderingin='"+orderingin+"'"
			if(self.verbose and jp.SysFrame(0,2)[-1][-2]==''): print
			return hp.ud_grade(inmap, nsideout, order_in=orderingin)





	def UnitConvert(self, unitin=None, unitout=None, inmap=None):
		''' '''
		unitin,unitout, k_in2out = self._UnitCheck(unitin,unitout)
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : strfont = ''
		else : strfont = '    '
		if (self.verbose) : 
			print strfont+'Tmap2Vis.UnitConvert'
			print strfont+strfont+"unitin='"+unitin+"', unitout='"+unitout+"'"
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : print
		if (inmap is None) : self.Tmap = self.Tmap * k_in2out
		else : return inmap * k_in2out





	def FreqinSet( self, freqin ) : 
		''' freqin must be in MHz '''
		self.freqin = float(freqin)
		if (self.verbose) : 
			print 'Tmap2Vis.FreqinSet'
			print '    freqin = %.3f MHz' % self.freqin



################### HealpixMap END #####################



	def InterArray( self, blorder, feedpos=None, telescope=None, vistype='cross1' ) : 
		'''
		feedpos:
			in meter
			feedpos = [(x1,y1,z1), (x2,y2,z2), ...]
			     North             +y
			West       East     -x    +x     zenith is +z
		        South             -y
		'''
		antarray = AntArray(telescope=telescope, verbose=False)
		antarray.verbose = self.verbose
		antarray.BlorderSet(blorder, feedpos)
		antarray.SelectVisType(vistype)
		self.baseline =antarray.Blorder.baseline[antarray.visorder]
		self.lon, self.lat = antarray.Ant.lonlat
		self.antform = antarray.Ant.antform





	def Pointing( self, pointDec ) : 
		'''
		pointDec:
			in degree
			Which Dec (Equational) does the antenna array points to?
			Can be one or N-D array (any shape)
		'''
		self.pointing = jp.npfmt(pointDec)
		if (self.verbose) : 
			print 'Tmap2Vis.Pointing'
			print '   ', self.pointing.round(3), 'degree \n'





	def Baseline2AntXYZ( self, baseline=None, lat=None ) : 
		'''
		baseline-xyz (local):
			     North             +y
			West       East     -x    +x     zenith is +z
		        South             -y
			
		AntXYZ:
			+Z: North Pole
			+X = +x
			+Y = +y

		baseline: in meter
		lat: in degree
		'''
		if (self.verbose) : print 'Tmap2Vis.Baseline2AntXYZ \n'
		coordtrans = jp.CoordTrans()
		if (baseline is None) : self.blXYZ = coordtrans.xyzRotation(self.baseline.T, ax=-(90-self.lat)).T
		else: return coordtrans.xyzRotation(baseline.T, ax=-(90-lat)).T




	def BeamMap( self, Dec, fwhm, nside=None ) : 
		'''
		Dec:
			in degree
			Dec of the antenna pointing
			Must be one, NOT ndarray

		fwhm:
			in degree
			FWHM of the beam
			Must be one, NOT ndarray

		return:
			Bmap is the beam of the antenna pointing to Dec in healpix
		'''
		if (nside is None) : nside = self.nside
		nside = jp.npfmt(nside, int).flatten()[0]
		Dec = jp.npfmt(Dec, float).flatten()[0]
		fwhm = jp.npfmt(fwhm, float).flatten()[0]
		freq = 1.03 *300/fwhm /self.antform *180/np.pi
		PrintStr = 'Dec='+('%.3f' % Dec)+' degree, fwhm='+('%.3f' % fwhm)+' degree, freq='+('%.3f' % freq)+' MHz'
		if (self.verbose) : 
			if (jp.SysFrame(0,2)[-1][-2] == '') : 
				print '    Tmap2Vis.BeamMap'
				print '        '+PrintStr
			else : print '    '+PrintStr
		Bmap = np.zeros(12*nside**2)
		npixc = hp.ang2pix(nside, (90-Dec)*np.pi/180, 0)
		Bmap[npixc] = 10000
		Bmap = hp.smoothing(Bmap, fwhm*np.pi/180, verbose=False)
	#	Bmap = Bmap / Bmap.max()
		Bmap = Bmap / Bmap.sum()
		Bmap[Bmap<0] = 0
		self.Bmap = Bmap  #@
		if (self.verbose and jp.SysFrame(0,2)[-1][-2]=='') : print
		return Bmap





	def Hpix2CeleXYZ( self, nside=None ) : 
		''' '''
		if (nside is None) : nside = self.nside
		else : nside = jp.npfmt(nside, int).flatten()[0]
		theta, phi = hp.pix2ang(nside, np.arange(12*nside**2))
		# RAhp, Dechp in Celestial-XYZ
		hpX = np.sin(theta) * np.cos(phi)
		hpY = np.sin(theta) * np.sin(phi)
		hpZ = np.cos(theta)
		self.hpixCeleXYZ = np.concatenate([hpX[:,None], hpY[:,None], hpZ[:,None]], -1)
		return self.hpixCeleXYZ





	def BeamMapComplex( self, RAnow, Bmap, freq ) : 
		''' return shape = (1, Nv, 12*nside**2) '''
		freq = jp.Num(freq, float)
		RAnow = jp.Num(RAnow, float)
		Bmap = jp.npfmt(Bmap)
		coordtrans = jp.CoordTrans()
		nside = hp.get_nside(Bmap)
		theta, phi = hp.pix2ang(nside, np.arange(12*nside**2))
		#--------------------------------------------------
		# For Bmap, rotate points in CeleXYZ from RA=0(+x) to RA=RAnow. Rotate points is negetive angle, so az=-RAnow
		az = -RAnow
		if (az != 0) : 
			pixnowB = hp.ang2pix(nside, theta, phi+az*np.pi/180)
			Bmap = Bmap[pixnowB]
			pixnowB = 0 #@
		# Now Bmapnow is still in CeleXYZ
		#--------------------------------------------------
		# For Tmap, don't need to rotate its points
		# Then for Tmap and Bmap(now), convert from CeleXYZ to AntXYZ. Note that here rotate XYZ, NOT points, so az=RAnow+90 from lon=0 to +AntX
		az = 90 + RAnow
		if (az != 0) : hpixXYZ = coordtrans.xyzRotation(self.hpixCeleXYZ.T, az=az).T
		else : hpixXYZ = self.hpixCeleXYZ.copy()
		#--------------------------------------------------
		# phase = 2pi/c * \vec(B) * \vec(s)
		phase = (2*np.pi/300*freq * self.blXYZ[:,None,:] * hpixXYZ[None,:,:]).sum(-1)  # (Nv, 12*nside**2) for 1 freq
		hpixXYZ = 0 #@
		#--------------------------------------------------
		Bmap = Bmap[None,:] * np.exp(1j*phase)
		return Bmap[None,:]  # None for Nfreq=1





	def Observate( self, freq, RAstart=None, obstimestart=None, inttime=None, Ntime=None, deltime=None, Nave=1 ) : 
		'''
		freq:
			in MHz
			Must be one, NOT N-D array, because when change freq, Bmap and Tmap will also change. Bmap: fwhm changes, Tmap: scale/K changes

		RAstart / obstimestart:
			Use RAstart first, if not, then use obstimestart
			If RAstart=obstimestart=None, set RAstart=0
			(1) RAstart:
				When start the observation, local meridian corresponds which RA?
				then get obstimestart
			(2) obstimestart:
				The start time of this observation in ephem format
				For example, '2016/11/05 03:17:32'
				then get RAstart

		inttime, Ntime, deltime: 
			set two of them

		inttime:
			in second
			Integration time of one output data point

		Ntime:
			deltime = Ntime * inttime

		deltime:
			in second
			Duration of this observation
			deltime = Ntime * inttime

		Nave:
			Every inttime output one data point.
			In this inttime, there are how many observated points averaged to one output data point.
		'''
		self.freq = jp.Num(freq, float)
		self.fwhm = 1.03 *300/freq /self.antform *180/np.pi
		if (inttime is not None and Ntime is not None) : 
			self.inttime = jp.Num(inttime, float)
			self.Ntime = jp.Num(Ntime, int)
			self.deltime = self.Ntime * self.inttime
		elif (inttime is not None and deltime is not None) : 
			self.inttime = jp.Num(inttime, float)
			self.deltime = jp.Num(deltime, float)
			self.Ntime = int(self.deltime / self.inttime)
			self.deltime = self.Ntime * self.inttime
		else : 
			self.Ntime = jp.Num(Ntime, int)
			self.deltime = jp.Num(deltime, float)
			self.inttime = self.deltime / self.Ntime
		if (Nave %2 != 1) : Nave += 1
		if (Nave < 1) : Nave = 1
		if (Nave == 1) : self.timelist = self.inttime/2 + self.inttime * np.arange(self.Ntime)
		else : self.timelist = np.linspace(0, self.deltime, (Nave-1)*self.Ntime+1)
		self.Nave = Nave
		timenow = jp.Time(0)
		RAnow = jp.SiderealTime(self.lon, date=timenow)*180/np.pi
		secnow = jp.Time(timenow)
		if (RAstart is None and obstimestart is None) : 
			self.RAstart = 0
			if (RAnow == 0) : sec = secnow
			else : sec = secnow + (360-RAnow)/15.*3600
			self.obstimestart = jp.Time(sec)
		elif (RAstart is not None) : 
			self.RAstart = jp.Num(RAstart, float)
			dRA = RAnow - RAstart
			if (dRA == 0) : sec = secnow
			elif (dRA > 0) : sec = secnow + (360-dRA)/15.*3600
			else : sec = secnow - dRA/15.*3600
			self.obstimestart = jp.Time(sec)
		else : 
			self.obstimestart = str(obstimestart)
			self.RAstart = jp.SiderealTime(self.lon, date=self.obstimestart)*180/np.pi
		if (self.verbose) : 
			print 'Tmap2Vis.Observate'
			print '    freq='+('%.3f' % self.freq)+' MHz, fwhm='+('%.3f' % self.fwhm)+' degree'
			print '    RAstart='+('%.3f' % self.RAstart)+' degree, obstimestart='+self.obstimestart
			print '    inttime='+str(self.inttime)+'sec, deltime='+('%.1f' % self.deltime)+' sec'
			print '    Nave='+str(self.Nave)+', Ntime='+str(self.Ntime)+'\n'





	def Vis( self, fwhm=None ) : 
		'''
		self.vis.shape = (pointing.size, Ntime, 1, bl.size)

		NOTE THAT here just generate 1 frequency, set in Tmap2Vis.Observate()
		'''
		if (self.verbose) : print 'Tmap2Vis.Vis'
		nside = self.nside
		coordtrans = jp.CoordTrans()
		if (fwhm is None) : freq, fwhm = self.freq, self.fwhm
		else : 
			fwhm = jp.npfmt(fwhm, float).flatten()
			freq = 1.03 *300/fwhm /self.antform *180/np.pi
		#--------------------------------------------------
		self.vis = []
		for i in xrange(self.pointing.size) : 
		#	verbose = self.verbose
			verbose = True
			Bmap = self.BeamMap(self.pointing[i], fwhm, nside)
			send = self.timelist
			bcast = (self.RAstart, Bmap, freq, self.hpixCeleXYZ, self.blXYZ, self.Tmap, verbose)
			#--------------------------------------------------
			if (self.Nprocess == 1) : 
				iterable = [(0,len(self.timelist)), send, bcast]
				vis = _DoMultiprocess_Vis(iterable)
			else : 
				pool =jp.PoolFor(0, len(self.timelist),self.Nprocess)
				vis =pool.map_async(_DoMultiprocess_Vis, send, bcast)
				vis = np.concatenate(vis, 0)
			if (self.verbose) : print
			#--------------------------------------------------
			if (self.Nave != 1) : 
				visa = []
				for j in xrange(self.Ntime) :  # Average
					n1 =    j *(self.Nave-1)
					n2 = (1+j)*(self.Nave-1) +1
					visa.append( vis[n1:n2].mean(0) )
				vis = np.array(visa)
				del visa
			#--------------------------------------------------
			self.vis.append( vis )
		self.vis = np.array(self.vis)
		if (self.verbose) : print



########################################################



	def SetLmax( self, lmax ) : 
		self.lmax = jp.Num(lmax, int)
		if (self.verbose) : 
			print 'Tmap2Vis.SetLmax'
			print '    lmax='+str(self.lmax)+'\n'



	def BeamLmax( self, verbose=None ) : 
		'''
		Decompose BeamMapComplex (with phase) to lm with lmax

		return:
			self.beamlmax.shape = (Npointing, Nv, (lmax+1)**2)
		'''
		if (self.verbose) : print 'Tmap2Vis.BeamLmax'
		if (verbose is None) : verbose = self.verbose
		lmax = self.lmax
		l = np.arange(lmax+1)
		# Negetive
		mn = np.arange(-l.max(), 0)
		idxdsn, morder = jp.LM2Index('almds', lmax, l,  mn)[:2]
		idxmn          = jp.LM2Index('alm'  , lmax, l, -mn)[0]
		# Positive
		mp = np.arange(0, l.max()+1)
		idxdsp = jp.LM2Index('almds', lmax, l, mp)[0]
		idxmp  = jp.LM2Index('alm'  , lmax, l, mp)[0]
		#--------------------------------------------------
		self.beamlmax = []
		for i in xrange(self.pointing.size) : 
			verbose0 = self.verbose
			self.verbose = verbose
			Bmap = self.BeamMap(self.pointing[i], self.fwhm)
			Bmap = self.BeamMapComplex(0, Bmap, self.freq)[0]
			self.verbose = verbose0
			beamlmtmp = []
			#--------------------------------------------------
			for j in xrange(Bmap.shape[0]) : 
				almR = hp.map2alm(Bmap[j].real, lmax)
				almI = hp.map2alm(Bmap[j].imag, lmax)
				almds = np.zeros((lmax+1)**2, complex)  # =beamlm
				#--------------------------------------------------
				almds[idxdsn] = (-1)**morder * (np.conj(almR[idxmn]) + 1j*np.conj(almI[idxmn]))
				almds[idxdsp] = almR[idxmp] + 1j*almI[idxmp]
				#--------------------------------------------------
				beamlmtmp.append( almds )
			self.beamlmax.append(np.array(beamlmtmp)) # (Nv, lmax)
		self.beamlmax = np.array(self.beamlmax)
		if (self.verbose) : print
		return self.beamlmax  # (Npoint, Nv, (lmax+1)**2)





	def BeamLM( self ) : 
		'''
		V(m)   = (-1)^m * Beam(l,-m) * Tmap(l,m)
		V*(-m) =          Beam*(l,m) * Tmap(l,m)
		(conj)              (conj)

		self.beamlm = [ (-1)^m *Beam(l,-m), Beam*(l,m) ]
		                                      (conj)
		self.beamlm.shape = (2, Npoint, Nv, m(lmax+1), l(lmax+1))
			2: [V(m), V*(-m)]

		NOTE THAT:
			Last two axes are m,l (l is the last one), NOT l,m
		'''
		if (self.verbose) : print 'Tmap2Vis.BeamLM \n'
		lmax = self.lmax
		l = np.arange(lmax+1)
		m = l.copy()
		#--------------------------------------------------
		# For Tmap ('Alm')
		valid  = True - jp.Invalid(jp.LM2Index('Alm', lmax, l, m, '2D')[0]).mask.T  # 2D, row-m, col-l
		#--------------------------------------------------
		# Vij(m), for beamlm ('AlmDS')
		idxBp = jp.LM2Index('AlmDS', lmax, l, -m, '2D')[0]
		idxBp = jp.Invalid(idxBp, 0).astype(int).T  # row-m, col-l
		# Vij^{*}(-m)
		idxBn = jp.LM2Index('AlmDS', lmax, l, m, '2D')[0]
		idxBn = jp.Invalid(idxBn, 0).astype(int).T  # row-m, col-l
		#--------------------------------------------------
		self.beamlm = np.array([(-1)**m *self.beamlmax[:,:,idxBp] *valid, np.conj(self.beamlmax[:,:,idxBn])*valid])
		return self.beamlm # (2, Npoint, Nv, m(lmax+1), l(lmax+1))





	def Tmaplm( self, what=None ) : 
		'''
		(1) what == None
			almT_1D = hp.map2alm(Tmap, lmax)
			return almT_1D

		(2) what is np.ndarray(almT_2D)
			Means what is almT_2D with shape=(l,m)
				NOTE THAT first axis(row) is l

			B^{-1} * V(m) = B^{-1} * B * almT_2D = almT_2D
			almT_2D = almT_1D[ almT_index ]
			return almT_1D

		You can obtain the reconstructed map by:
			Tmap_re = hp.alm2map(self.Tmaplm, nside)
		'''
		if (self.verbose) : 
			print 'Tmap2Vis.TmapLM \n'
			if (what is None) : print '    what=None'
			else : print '    what.shape='+str(what.shape)+' array'
			print
		#--------------------------------------------------
		if (which is None) : 
			self.Tmaplm = hp.map2alm(self.Tmap, self.lmax)
			return self.Tmaplm
		#--------------------------------------------------
		lmax = self.lmax
		nside = lmax/4.
		n = 1+int(np.log(nside) / np.log(2))
		nside = 2**n
		self.Tmaplm = hp.map2alm(np.zeros(12*nside**2), lmax)
		#--------------------------------------------------
		l = np.arange(lmax+1)
		m = l.copy()
		#--------------------------------------------------
		# For Tmap ('Alm')
		index = jp.LM2Index('Alm', lmax, l, m, '2D')[0]
		valid = True - jp.Invalid(index).mask
		index = index[valid].astype(int)  # 1D
		self.Tmaplm[index] = what[valid]
		return self.Tmaplm





	def VisM( self ) : 
		'''
		Use BeamLM and TmapLM to calculate V(m)/vism
		vism is the FFT of vis

		return:
			self.vism = [V(m), V(-m)]  # NOT conj
			self.vism.shape = (2, Npoint, 1, Nv, lmax+1)
		                                 Nfreq=1
		'''
		if (self.verbose) : print 'Tmap2Vis.VisM \n'
		lmax, shape = self.lmax, self.beamlm.shape
		Np, Nv, Nl = shape[1], shape[2], shape[4]
		self.vism = np.zeros([2, Np, 1, Nv, Nl], complex)
		almT = hp.map2alm(self.Tmap, lmax)
		#--------------------------------------------------
		l = np.arange(lmax+1)
		m = l.copy()
		# For Tmap ('Alm')
		idxT = jp.LM2Index('Alm', lmax, l, m, '2D')[0]
		idxT = jp.Invalid(idxT , 0).astype(int).T  # ml
		#--------------------------------------------------
		self.vism[0,:,0,:,:] = (self.beamlm[0] * almT[idxT][None,None,:,:]).sum(-1)
		self.vism[1,:,0,:,:] = np.conj( (self.beamlm[1] * almT[idxT][None,None,:,:]).sum(-1) )
		return self.vism





	def VisM2Vis( self ) : 
		''' 
		Use self.vism to obtain self.vis2 by ifft
		self.vism2.shape = (Npoint, 2*lmax+1, 1, Nv)
		'''
		vm, v2 = self.vism
		v2 = v2[:,:,:,:-1][:,:,:,::-1]
		vm = np.concatenate([vm, v2], -1)
		del v2
		self.vis2 = np.fft.ifft(vm)
		del vm
		self.vis2 = jp.ArrayAxis(self.vis2, -1, 1, 'move')
		return self.vis2




