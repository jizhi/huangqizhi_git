import h5py
import pyfits
import sys
import ephem
import numpy as np
import jizhipy as jp
from jizhipy.Plot import *
##################################################



class _AntArrayBlorder( object ) : 
	''' self.Blorder '''


	def _Blorder( self, blorder ) : 
		try : self.blorder = blorder[:]
		except : jp.Raise(Exception, 'Missing AntArray.Blorder.blorder')
		self._blorder = self.blorder.copy()



	def _Feedpos( self, feedpos ) : 
		self.feedpos = feedpos[:]
		shape = self.feedpos.shape
		self.channelpos = np.concatenate([self.feedpos[:,None,:], self.feedpos[:,None,:]], 1).reshape(2*shape[0], shape[1])



	def _AutoCross( self ) : 
		''' Return the index of auto, cross with H-H and V-V, cross with H-V '''
		blorder = self.blorder
		dblorder = blorder[:,0] - blorder[:,1]  # For auto
		nblorder = np.arange(len(blorder))
		#----------
		tfa = (dblorder==0)
		nauto = nblorder[tfa]
		blauto = blorder[tfa,0].copy() %2
		self.auto1 = nauto[blauto==1] # H
		self.auto2 = nauto[blauto==0] # V
		#----------
		dblorder = abs(blorder[:,0]%2 - blorder[:,1]%2)
		tfc1 = (dblorder==0) * (True-tfa)
		ncross12 = nblorder[tfc1]
		blcross12 = blorder[tfc1,0].copy() %2
		self.cross1 = ncross12[blcross12==1] # H
		self.cross2 = ncross12[blcross12==0] # V
		#----------
		tfc2 = (dblorder==1) * (True-tfa)
		self.cross3 = nblorder[tfc2]



	def _Baseline( self ) : 
		blorder = self.blorder.copy()
		feedpos = self.feedpos.copy()
		blorder[blorder%2==1] = blorder[blorder%2==1]+1
		blorder = blorder/2-1
		feedpos1, feedpos2 = feedpos[blorder[:,0]], feedpos[blorder[:,1]]
		self.baseline = feedpos2 - feedpos1



	def Bl2Order( self, *arg ) : 
		''' Bl2Order(1,3) or Bl2Order((1,3)), return order= '''
		if (len(arg) == 1) : bl = arg[0]
		else : bl = arg[:2]
		bl = np.sort(bl)
		dblorder = abs(self.blorder-np.array(bl)[None,:]).sum(1)
		order = np.arange(dblorder.size)[dblorder==0]
		if (order.size == 0) : Raise(Exception, 'NOT exisit baseline='+str(bl[0])+'-'+str(bl[1]))
		return order[0]


	def Order2Bl( self, order ) : 
		return self.blorder[int(round(order))]



##################################################
##################################################
##################################################



class _AntArrayFile( object ) : 
	''' self.File '''



	def _FileList( self, filedir, telescope, filefmt ) : 
		''' self.filedir, self.filelist '''
		self.filedir = jp.AbsPath(filedir)
		self.filedirname = self.filedir.split('/')[-2]
		if (telescope == 'tianlai') : self.filefmt = '.hdf5'
		elif (telescope == 'paon4') : self.filefmt = '.fits'
		else : self.filefmt = '.'+str(filefmt).split('.')[-1]
		if (self.filefmt in ['.None', '.']) : jp.Raise(Exception, "filefmt in [None, '']")
		old = '_old'+self.filefmt
		if (self.filefmt == '.hdf5') : 
			filelist = jp.ShellCmd('ls '+self.filedir+'*.hdf5')
		elif (self.filefmt == '.fits') : 
			filelist = jp.ShellCmd('ls '+self.filedir+'*.fits')
		ndel = 0
		for i in xrange(len(filelist)) : 
			if (filelist[i-ndel][-len(old):] == old) : 
				filelist.pop(i-ndel)
				ndel += 1
		self.filelist = filelist
		if (len(filelist) == 0) : print '    Warning: self.filelist == [], is empty'



	def _WhichFile( self, which, telescope ) : 
		# self.File.nfile
		which = np.sort(jp.npfmt(which).round().astype(int).flatten())
		self.nfile = which[(0<=which)*(which<len(self.filelist))]
		which = np.arange(self.nfile[0], self.nfile[-1]+1)
		if (self.nfile.size != which.size) : 
			missing = []
			for i in xrange(which.size) : 
				if (which[i] not in self.nfile) : missing.append(which[i])
			print '    Warning: self.File.nfile is from '+str(which[0])+'(include) to '+str(which[-1])+'(include), but now missing '+str(missing)
		#--------------------------------------------------
		# self.File.filepath
		self.filepath = list(np.array(self.filelist)[self.nfile])
		#--------------------------------------------------
		if (telescope == 'tianlai') : 
			fo = h5py.File(self.filepath[0], 'r')
			self.obstime = fo.attrs['obstime']
			self.sec1970 = fo.attrs['sec1970']
		elif (telescope == 'paon4') : 
			fo = pyfits.open(self.filepath[0])
			obstime = fo[0].header['dateobs']
			for i in xrange(len(obstime)) : 
				if (obstime[i] == '-') : obstime = obstime[:i]+'/'+obstime[i+1:]
				elif (obstime[i] == 'T') : obstime = obstime[:i]+' '+obstime[i+1:]
			self.obstime = obstime
			self.sec1970 = jp.Time(self.obstime, 0)



##################################################
##################################################
##################################################



class _AntArrayAnt( object ) : 
	''' self.Ant '''
	pass



##################################################
##################################################
##################################################



class AntArray( object ) : 
	''' 
	self.MaskChannel(), self.SelectChannel():
		Modify self.Blorder.maskorder, maskchannel, selectorder, selectchannel, auto1, auto2, cross1, cross2, cross3

	self.SelectVisType():
		Modifies self.visorder

	self.Blorder.blorder, baseline, channelpos won't change at any time !!!

	maskorder, selectorder, auto1, auto2, cross1, cross2, cross3, visorder

	maskchannel, selectchannel:
		Are used to channelpos[selectchannel]
	'''



	def __init__( self, filedir=None, telescope=None, filefmt=None, Nprocess=None, verbose=True, outdir='' ) :
		'''
		filedir:
			In this directory/folder (absolute path), all .hdf5/.fits are one continuum observation splited into several files. If one file one observation, this filedir must just contain one file.

		telescope:
			'tianlai' | 'paon4' | None

		Nprocess:
			number of multiprocess

		verbose:
			True | False
			Whether print the messages

		##################################################

		blorder:
			blorder with acORall='all'
			For PAON-4, len(blorder) = 36

		self.__dict__.keys() 
			= [WhichFile(), MaskChannel(), SelectChannel(), SelectVisType(), Blorder, File, Ant, vis, vistype, visorder]

		self.Blorder.__dict__.keys() 
			blorder starts from 1, not 0 !
			= [blorder, auto1, auto2, cross1, cross2, cross3, feedpos, channelpos, baseline, Bl2Order(), Order2Bl(), maskchannel, maskorder, selectchannel, selectorder]
	
		self.File.__dict__.keys() 
			= [filedir, filelist, filepath, nfile, obstime, sec1970,  nfiletransit, transittimelocal, transittimetotal]
				* transittimelocal, transittimetotal are pixels, NOT time in second

		self.Ant.__dict__.keys()
			= [noisesourcepos, inttime, lonlat, freq, antform]

		transittime:
			in second, start from filelist[0] (the beginning)
			If NO fo['transitsource'], set self.File.nfiletransit and self.File.transittime (second from filelist[0]) here
		'''
		self.starttime = jp.Time(1)
		self.verbose, self.Nprocess = bool(verbose), Nprocess
		if (self.verbose) : print '-------------------- AntArray --------------------\n'
		self.Nprocess = jp.NprocessCPU(Nprocess, verbose)[0]
		if (telescope is not None) : 
			if (str(telescope).lower() not in ['tianlai', 'paon4']) : jp.Raise(Exception, 'telescope="'+str(telescope)+'" not in ["tianlai", "paon4"]')
		telescope = str(telescope).lower()
		self.Blorder = _AntArrayBlorder()
		self.File    = _AntArrayFile()
		self.Ant     = _AntArrayAnt()
		self.Ant._telescope = telescope
		#--------------------------------------------------
		if (telescope == 'paon4') : 
			self.Ant.lonlat = np.array([2.1996, 47.3820])
			self.Ant.telescope = 'PAON-4'
			self.Ant.antform = 5.0
			self.Ant.timezone = 0  # obstime in PAON4 is UTC
		if (filedir is None) : return
		#--------------------------------------------------
		self.File._FileList(filedir, telescope, filefmt)
		if (self.File.filefmt == '.hdf5') : 
			fo = h5py.File(self.File.filelist[0], 'r')
		elif (self.File.filefmt == '.fits') : 
			fo = pyfits.open(self.File.filelist[0])
		#--------------------------------------------------
		self._Ant(fo)
		self._File(fo)
		self._Blorder(fo, None)
		verbose = self.verbose
		try : 
			self.verbose = False
			self.MaskChannel()
			self.SelectVisType('all')
		except : pass
		self.verbose = verbose
		self.outdir = jp.Outdir((None,'file'), (0,'file'))
		try : self.outdir = jp.Mkdir(self.outdir + outdir)
		except : pass
		self.Ant.freqorder = None
		self.reducetime = 1
		self.reducefreq = 1
		self.outname = self.outdir + 'AntArray.vis_'+self.File.filedirname+'.hdf5'





	def _Ant( self, fo ) : 
		if (self.Ant._telescope == 'tianlai') : 
			self.Ant.lonlat = np.array([fo.attrs['sitelon'], fo.attrs['sitelat']])
			self.Ant.inttime = fo.attrs['inttime']*1.
			self.Ant.inttimereal = self.Ant.inttime
			self.Ant.telescope = fo.attrs['telescope']
			try : self.Ant.noisesourcepos = fo['noisesource'][0]
			except : pass
			self.Ant.siteelev = fo.attrs['siteelev']
			self.Ant.timezone = int(fo.attrs['timezone'][3:-1])
			self.Ant.N0 = fo['vis'].shape[0]
			#--------------------------------------------------
			if ('cylinder' in self.Ant.telescope.lower()) : 
				self.Ant.antform = np.array([fo.attrs['cywid'], fo.attrs['cylen']])
			else : self.Ant.antform = fo.attrs['dishdiam']
			#--------------------------------------------------
			freq1 = fo.attrs['freqstart'] + np.arange(fo.attrs['nfreq']) * fo.attrs['freqstep']
			freq2 = np.linspace(fo.attrs['freqstart'], fo.attrs['freqstart']+fo.attrs['freqstep']*fo.attrs['nfreq'], fo.attrs['nfreq'])
			self.Ant.freq = freq1 if(freq1.size==fo.attrs['nfreq'])else freq2
			#--------------------------------------------------

		elif (self.Ant._telescope == 'paon4') : 
			hdr = fo[0].header
			self.Ant.N0 = fo[0].header['naxis3']
			self.Ant.inttime = 1.*hdr['deltime'] / hdr['naxis3']
			self.Ant.inttimereal = hdr['npaqsum'] * 8192/500e6 #@
			self.Ant.freq = np.linspace(1250, 1500, hdr['naxis1'])





	def _File( self, fo ) : 
		'''
		self.File.transittime
			in second, NOT pixel, start from the begining of this observation (self.File.obstimestart)
		'''
		if (self.Ant._telescope == 'tianlai') : 
			self.File.obstimestart = fo.attrs['obstime']
			self.File.sec1970start = fo.attrs['sec1970']
			try : # transitsource
				tsa = np.array(fo['transitsource'][:], float)
				self.File.transitname = fo['transitsource'].attrs.items()[-1][-1].split(',')
				self.File.transitRADec = tsa[:,1:3]
				self.File.transitAzAlt = tsa[:,3:]
				transittime = tsa.T[0]
				if (transittime.size==1 and transittime[0]<0) : transittime[0] = 0
				elif (transittime[transittime<0].size > 0) : 
					transittime[transittime<0] = transittime[transittime>=0][0]+np.zeros(transittime[transittime<0].size)
				# Recal RADec if is Sun
				for i in xrange(len(self.transitname)) : 
					name = self.transitname[i].lower()
					if (name in ['sun', 'solar']) :
						gatech = ephem.Observer()
						gatech.lon = str(self.Ant.lonlat[0])
						gatech.lat = str(self.Ant.lonlat[1])
						gatech.elevation = self.Ant.siteelev
						gatech.date = jp.Time(jp.Time(transittime[0], self.Ant.timezone), self.timezone, 0)
						sun = ephem.Sun(gatech)
						self.File.transitRADec[i] = [sun.a_ra*180/np.pi, sun.a_dec*180/np.pi]
				self.File.transittime =transittime-self.sec1970start
				self.File.nfiletransit = (self.transittime / Ant.inttime / Ant.N0).astype(int)
			except : print "    Warning: NOT exist fo['transitsource'], you can set it by using AntArray.TransitsourceSet()\n"
			#--------------------------------------------------

		elif (self.Ant._telescope == 'paon4') : 
			obstime = fo[0].header['dateobs']  # UTC
			for i in xrange(len(obstime)) : 
				if (obstime[i] == '-') : obstime = obstime[:i]+'/'+obstime[i+1:]
				elif (obstime[i] == 'T') : obstime = obstime[:i]+' '+obstime[i+1:]
			self.File.obstimestart = obstime
			self.File.sec1970start = jp.Time(self.File.obstimestart, 0)
			print "    Warning: NOT exist fo['transitsource'], you can set it by using AntArray.TransitsourceSet()\n"





	def TransitsourceSet( self, transitsource ) : 
		'''
		transitsource[i] is for i-th source
		transitsource[i] = [
			transitname,  (str, like 'Sun')
			transittime,  (second, NOT pixel, starts from the starting time of this observation)
			RA,           (degree)
			Dec]          (degree)
		For example:
			transitsource = [('CygA', 13540.6, 299.8683, 40.7339)]
		'''
		if (self.verbose) : print 'AntArray.TransitsourceSet'
		istype = jp.IsType()
		if (istype.isstr(transitsource[0])) : transitsource = [transitsource]
		self.File.transitname = []
		for i in xrange(len(transitsource)) : 
			if (self.verbose) : print '    ('+transitsource[i][0]+', %.3f, %.3f, %.3f)' % tuple(transitsource[i][1:])
			self.File.transitname.append(str(transitsource[i][0]))
			transitsource[i] = transitsource[i][1:]
		transitsource = jp.npfmt(transitsource)
		self.File.transittime = transitsource[:,0]
		self.File.transitRADec = transitsource[:,1:]
		self.File.nfiletransit = (self.File.transittime / self.Ant.inttime / self.Ant.N0).astype(int)
		if (self.verbose) : print





	def _Blorder( self, fo, blorder, feedpos=None ) : 
		if (blorder is None) : 
			if (self.Ant._telescope == 'tianlai') : 
				blorder = fo['blorder']
			elif (self.Ant._telescope == 'paon4') : 
				print "    Warning: NOT exist fo['blorder'], you can set it by using AntArray.BlorderSet()\n"
				return
		else : blorder = jp.npfmt(blorder)
		if (feedpos is None) : 
			if (self.Ant._telescope == 'tianlai') : 
				feedpos = fo['feedpos']
			elif (self.Ant._telescope == 'paon4') : 
				# PAON-4
				#			4         North            +y
				#	2	1        West     East     -x    +x
				#			3         South            -y
				feedpos = np.array([[0,0,0], [-5.993,-0.001,0], [4.38,-5.996,0], [4.383,5.995,0]])
		else : feedpos = jp.npfmt(feedpos)
		self.Blorder._Feedpos(feedpos)
		self.Blorder._Blorder(blorder)
		self.Blorder._AutoCross()
		self.Blorder._Baseline()





	def BlorderSet( self, blorder, feedpos=None ) : 
		'''
		PAON-4
					4         North            +y
			2	1        West     East     -x    +x
					3         South            -y
		'''
		if (self.Ant._telescope == 'tianlai') : return
		if (self.verbose) : 
			print 'AntArray.BlorderSet'
			printstr = ''
			for i in xrange(len(blorder)) : printstr += '(%i,%i),' % tuple(blorder[i])
			print '    ['+printstr[:-1]+']'
		self._Blorder(None, blorder, feedpos)
		verbose, self.verbose = self.verbose, False
		self.MaskChannel()
		self.SelectVisType('all')
		self.verbose = verbose
		if (self.verbose) : print





	def WhichFile( self, whichfile ) : 
		'''
		whichfile:
			(1) Can be absolute path: whichfile = '/xxx/yyy.hdf5'
			(2) Or index/number of the file: whichfile = 3
		**	(1,2) whichfile can be one or list: 
				whichfile = [0,2,5,7,9,13,18]
				whichfile = ['a.hdf5', 'k.hdf5', 'x.hdf5']
				If so, will merge them into one array one by one!
			(3) string with 'transitsource':
				1. =='transitsource': first transitsource
				2. =='transitsource5' or 'transitsource-5' or 'transitsource_5': 5-th transitsource (start from 1, not 0)
			(4) whichfile='all'
		'''
		if (self.verbose) : print 'AntArray.WhichFile'
		istype = jp.IsType()
		if ('transitsource' in str(whichfile).lower()) : 
			if ('nfiletransit' not in self.File.__dict__.keys()) : jp.Raise(Exception, "NOT exist fo['transitsource'], you can set it by using AntArray.Transitsource()")
			else : 
				n =str(whichfile[13:]).split('-')[-1].split('_')[-1]
				try : n = int(n)
				except : n = 0
				which = jp.npfmt(self.File.nfiletransit[n])
		elif (str(whichfile).lower() == 'all') : 
			which = np.arange(len(self.File.filelist))
		else : which = whichfile
		self.File._WhichFile(which, self.Ant._telescope)
		if (self.verbose) : 
			print '    whichfile='+str(list(self.File.nfile))+'\n'





	def MaskChannel( self, maskchannel=None ) : 
		'''
		self.maskorder, self.selectorder: for the whole 528 baselines, for self.vis[None,None,:]
		self.maskchannel, self.selectchannel: for the whole 32 channels

		maskchannel:
			Start from 1 (not 0)
			int or [] or () or np.array([])
			I use SunCygnusA_20160603 and CasA_20160702 to check the channels, very bad channels: [9,10,21,23,24,27], a little bad channel: [28], may be bad/good: [11]
		'''
		#--------------------------------------------------
		if (maskchannel is None) : maskchannel = []
		self.Blorder.maskchannel = np.sort(jp.npfmt(maskchannel))
		#--------------------------------------------------
		if (self.Blorder.maskchannel.size == 0) : 
			self.Blorder.maskorder = np.array([])
			self.Blorder.selectchannel = np.arange(len(self.Blorder.channelpos))+1
			self.Blorder.selectorder = np.arange(len(self.Blorder._blorder))
			return
		#--------------------------------------------------
		_maskorder = np.zeros(len(self.Blorder._blorder), bool)
		for i in self.Blorder.maskchannel : 
			_maskorder[self.Blorder._blorder[:,0]==i] = True
			_maskorder[self.Blorder._blorder[:,1]==i] = True
		self.Blorder.maskorder = np.arange(_maskorder.size)[_maskorder]
		#--------------------------------------------------
		_selectchannel = np.ones(len(self.Blorder.feedpos)*2, bool)
		_selectchannel[self.Blorder.maskchannel-1] = False
		self.Blorder.selectchannel = np.arange(len(self.Blorder.feedpos)*2)[_selectchannel]+1
		#--------------------------------------------------
		_selectorder = np.ones(len(self.Blorder._blorder),bool)
		_selectorder[self.Blorder.maskorder] = False
		self.Blorder.selectorder = np.arange(_selectorder.size)[_selectorder]
		#--------------------------------------------------
		self._MaskChannel()  # Recal auto, cross1, cross2
		verbose, self.verbose = self.verbose, False
		self.SelectVisType(self.vistype)  # Reset visorder
		self.verbose = verbose





	def SelectChannel( self, selectchannel=None ) : 
		''' selectchannel: Start from 1 (not 0) '''
		if (selectchannel is None) : selectchannel = []
		selectchannel = jp.npfmt(selectchannel)-1
		_maskchannel = np.ones(len(self.Blorder.feedpos)*2, bool)
		_maskchannel[selectchannel] = False
		maskchannel = np.arange(len(self.Blorder.feedpos)*2)[_maskchannel]+1
		self.MaskChannel(maskchannel)





	def _MaskChannel( self ) : 
		self.Blorder.blorder = self.Blorder._blorder[self.Blorder.selectorder]
		self.Blorder._AutoCross()
		self.Blorder.blorder = self.Blorder._blorder
		norder = np.arange(len(self.Blorder._blorder))
		self.Blorder.auto1 = norder[self.Blorder.selectorder][self.Blorder.auto1]
		self.Blorder.auto2 = norder[self.Blorder.selectorder][self.Blorder.auto2]
		self.Blorder.cross1 = norder[self.Blorder.selectorder][self.Blorder.cross1]
		self.Blorder.cross2 = norder[self.Blorder.selectorder][self.Blorder.cross2]
		self.Blorder.cross3 = norder[self.Blorder.selectorder][self.Blorder.cross3]





	def SelectVisType( self, vistype='cross1' ) : 
		self.vistype = str(vistype).lower()
		if (self.verbose) : 
			print 'AntArray.SelectVisType'
			print '    vistype='+self.vistype+'\n'
		if   (self.vistype == 'auto1') : 
			self.visorder = self.Blorder.auto1
		elif (self.vistype == 'auto2') : 
			self.visorder = self.Blorder.auto2
		elif (self.vistype == 'cross1') : 
			self.visorder = self.Blorder.cross1
		elif (self.vistype == 'cross2') : 
			self.visorder = self.Blorder.cross2
		elif (self.vistype == 'cross3') : 
			self.visorder = self.Blorder.cross3
		elif (self.vistype == 'all') : 
			self.visorder = self.Blorder.selectorder
		else : 
			if (self.verbose) : print "    Warning: self.vistype='"+vistype.lower()+"' not in ['auto1', 'auto2', 'cross1', 'cross2, cross3'], reset to self.vistype='cross1'\n"
			self.vistype = 'cross1'
			self.visorder = self.Blorder.cross1
		if (self.visorder.size == 0) : jp.Raise(Exception, 'self.visorder == [], NO data')





	def SelectFreq( self, nfreq, unit=None ) : 
		'''
		nfreq:
			int/float/list

		unit:
			None, 'Hz'/'MHz'/'GHz'
			==None: nfreq is number
			=='Hz'/'MHz'/'GHz': nfreq is frequency value

		nfreq, unit:
			If ('hz' not in str(unit).lower()) : 
					self.freqorder = nfreq
					antarray.vis[:, self.freqorder, self.visorder]
			Else, use nfreq to get self.freqorder
		'''
		nfreq = jp.npfmt(nfreq)
		if ('hz' not in str(unit).lower()) : 
			nfreq = nfreq.astype(int)
			unit = ''
		else : unit = unit.lower()
		#--------------------------------------------------
		if (self.verbose) : print 'AntArray.SelectFreq'
		if ('vis' in self.__dict__.keys()) : 
			jp.Raise(Warning, 'You should do AntArray.SelectFreq() before AntArray.Vis()\n         Do nothing !')
			return
		if (self.Ant.freqorder is not None) : 
			jp.Raise(Warning, 'You already did AntArray.SelectFreq() once OR did AntArray.ReduceVis(reducefreq != None)\n         Do nothing !')
			return
		#--------------------------------------------------
		if (self.verbose) : 
			if (unit == '') : 
				if (nfreq.size <= 12) : pstr = str(list(nfreq)) 
				else : 
					s1 = str(list(nfreq[:5]))[:-1]
					s2 = str(list(nfreq[-5:]))[1:]
					pstr = s1+', ..., '+s2
				print '    nfreq='+pstr+'\n'
			else : 
				pstr = '['
				if (nfreq.size <= 12) :
					for i in xrange(nfreq.size) : 
						pstr += ('%.3f, ' % nfreq[i])
					pstr = pstr[:-2]+']'
				else : 
					pstr1, pstr2 = '', ''
					f1, f2 = nfreq[:5], nfreq[-5:]
					for i in xrange(f1.size) : 
						pstr1 += ('%.3f, ' % f1[i])
						pstr2 += ('%.3f, ' % f2[i])
					pstr += pstr1 + '..., ' + pstr2[:-2] + ']'
				print '    freq='+pstr+' MHz \n'
		if (unit == '') : 
			self.Ant.freqorder = nfreq
			self.freq = None
		else : 
			if (unit == 'hz') : nfreq = nfreq * 1e-6  # MHz
			elif (unit == 'ghz') : nfreq = nfreq * 1e3
			else : nfreq = 1.*nfreq
			dfreq = abs(nfreq[:,None] - self.Ant.freq[None,:])
			self.Ant.freqorder = np.where(dfreq==dfreq.min(1)[:,None])[1]
			self.freq = self.Ant.freq[self.Ant.freqorder]





	def ReduceVis( self, reducetime=None, reducefreq=None ) : 
		''' 
		Reduce the vis matrix by averaging the bins

		self.freq:
			Same shape as self.vis.shape[1]

		(1) If AntArray.SelectFreq() with unit=None/''
			use existing self.Ant.freqorder and jp.Smooth(self.Ant.freqorder, reducefreq)

		(2) If AntArray.SelectFreq() with 'Hz' in unit
			recalculate freqorder and reduce it
		'''
		if (reducetime is not None and int(round(reducetime))>1) : 
			self.reducetime = int(round(reducetime))
		if (reducefreq is not None and int(round(reducefreq))>1) : 
			self.reducefreq = int(round(reducefreq))
			if (self.Ant.freqorder is None) : 
				self.Ant.freqorder = np.arange(self.Ant.freq.size) 
				self.freq = jp.Smooth(self.Ant.freq, 0, self.reducefreq, reduceshape=True)
			else : 
				if (self.freq is not None) : 
					freqorder = []
					sfo, N, Ntot = self.Ant.freqorder, self.reducefreq, self.Ant.freq.size
					for i in xrange(sfo.size) : 
						fo = np.arange(sfo[i]-N/2, sfo[i]+N/2+N%2)
						fo[fo<0] = 0
						fo[fo>=Ntot] = Ntot-1
						freqorder.append(fo)
					self.Ant.freqorder = np.concatenate(freqorder)
				self.freq = self.Ant.freq[self.Ant.freqorder]
				self.freq = jp.Smooth(self.freq, 0, self.reducefreq, reduceshape=True)
		if (self.verbose) : 
			print 'AntArray.ReduceVis'
			print '    reducetime='+str(self.reducetime)
			print '    reductfreq='+str(self.reducefreq)+'\n'





	def Vis( self, inhdf5=None, xhour=False ) : 
		'''
		Read the visibibity data

		inhdf5:
			(1) ==None, read data from self.File.filepath
			(2) path of 'AntArray.Vis_xxx.hdf5', read vis, freq, timem, bli, blo from this file

		self.vis:
			self.vis.shape = (time, freq, baseline)

		self.timem:
			in minute
			self.timem[0] is the time(minute) on day self.File.obstime.split(' ')[0] (start from 0:0:0 A.M.)
			Note that it is time at self.Ant.timezone!

		self.freq, self.blo, self.bli
		'''
		if (self.verbose) : print 'AntArray.Vis'
		if (inhdf5 is not None) : 
			if (self.verbose) : print "    Reading '"+inhdf5+"' \n"
			fo = h5py.File(inhdf5, 'r')
			self.vis = fo['vis'][:]
			self.freq = fo['freq'][:]
			self.blo = fo['blo'][:]
			self.bli = fo['bli'][:]
			self.File.transitRADec = fo['RADec'][:]
			self.Ant.lonlat = fo['lonlat'][:]
			self.timem = fo['timem'][:]
			self.Ant.timezone = fo['timem'].attrs['timezone']
			t = self.timem[0] / 60.
			h = int(t)
			m = int((t - h)*60)
			s = (t - h - m/60.)*3600
			h = str(h) if(h>=10)else '0'+str(h)
			m = str(m) if(m>=10)else '0'+str(m)
			s = ('0%.2f' % s) if(s>=10)else ('0%.2f' % s)
			self.File.obstime = fo['timem'].attrs['obsdate']+' '+h+':'+m+':'+s
			self.File.filedirname = fo.attrs['file'][:-1]
			self.vistype = fo.attrs['vistype']
			self.Ant._telescope = fo.attrs['telescope']
			fo.close()
		#--------------------------------------------------

		else : 
			if (self.verbose) : print '    start @', jp.Time(1)
			if (self.Ant.freqorder is None) : 
				self.Ant.freqorder = np.arange(self.Ant.freq.size) 
				self.freq = self.Ant.freq
			vistype = self.vistype[:-1]
			self.vis = []
			for i in xrange(len(self.File.filepath)) : 
				if (self.verbose) : print '    Reading file '+str(self.File.nfile[i])+' = "'+self.File.filepath[i].split('/')[-1]+'"'
				if (self.Ant._telescope == 'tianlai') : 
					fo = h5py.File(self.File.filepath[i], 'r')
					data = jp.Getdata(fo['vis'], None, self.Ant.freqorder, self.visorder)
					if (vistype =='auto') : data = data.real
					#----------------------------------------
	
				elif (self.Ant._telescope == 'paon4') : 
					fo = pyfits.open(self.File.filepath[i])
					if (len(fo.info(0)) == 2) : 
						data = jp.Getdata(fo[0].data, None, self.visorder, self.Ant.freqorder)
						if (vistype == 'cross'): 
							data = data + 1j*jp.Getdata(fo[1].data, None, self.visorder, self.Ant.freqorder)
					#----------------------------------------
	
					else : 
						dbl = self.Blorder.blorder[:,0] - self.Blorder.blorder[:,1]
						if (vistype == 'auto') : 
							order  = self.Blorder.blorder[dbl==0]
						elif(vistype == 'cross') : 
							order = self.Blorder.blorder[dbl!=0]
							tf = abs(order[:,0] - order[:,1]) % 2
							order = order[tf==0]
						# visorder of len(fo.info(0))==3
						self._visorder3 = []
						bl = self.Blorder.blorder[self.visorder]
						for i in xrange(len(bl)) : 
							try : self._visorder3.append( np.where((order[:,0]==bl[i,0])*(order[:,1]==bl[i,1]))[0][0] )
							except : pass
						self._visorder3 = np.array(self._visorder3)
						if (self._visorder3.size == 0) : jp.Raise(Exception, 'self._visorder3 == [], NO data')
						if (vistype == 'auto') : 
							data = jp.Getdata(fo[0].data, None, self._visorder3, self.Ant.freqorder)
						elif (vistype == 'cross'): 
							data = jp.Getdata(fo[1].data, None, self._visorder3, self.Ant.freqorder)
							data = data + 1j*jp.Getdata(fo[2].data, None, self._visorder3, self.Ant.freqorder)
				#--------------------------------------------------
	
				if (self.Ant._telescope == 'paon4') : axis = 2
				elif (self.Ant._telescope == 'tianlai') : axis = 1
				data = jp.Smooth(data, axis, self.reducefreq, reduceshape=True)
				self.vis.append( data )
				#--------------------------------------------------
			self.vis = np.concatenate(self.vis, 0)
			if (self.Ant._telescope == 'paon4') : 
				Nf, Nv = self.Ant.freqorder.size, self.visorder.size
				if (Nv <= Nf) : self.vis = jp.ArrayAxis(self.vis, 1, 2, 'move')
				else : self.vis =jp.ArrayAxis(self.vis, 2,1, 'move')
			#--------------------------------------------------
	
			localtime = np.array(self.File.obstime.split(' ')[-1].split(':'), float)
			localtime = localtime[0]*60+localtime[1]+localtime[2]/60
			self.timem = np.arange(len(self.vis)) * self.Ant.inttime/60 + localtime # real time in minute starting from 00:00
			#--------------------------------------------------
	
			self.vis = jp.Smooth(self.vis, 0, self.reducetime, reduceshape=True)
			self.timem = jp.Smooth(self.timem, 0, self.reducetime, reduceshape=True)
			self.blo = self.Blorder.blorder[self.visorder]
			self.bli = self.Blorder.baseline[self.visorder]
			if (self.verbose) : 
				print '\n    AntArray.vis.shape='+str(self.vis.shape)
#			self.PlotVis(None, None, xhour)
			if (self.verbose) : print '    end   @', jp.Time(1)+'\n'





	def Save( self, outname=None ) : 
		'''
		outname:
			Absolute path of output .hdf5
		'''
		if (outname is not None) : self.outname = str(outname)
		if (self.outname[-5:] !='.hdf5') : self.outname += '.hdf5'
		if (self.verbose) : print 'AntArray.Save'
		if (self.verbose) : print '    Saving to  "'+jp.AbsPath(self.outname)+'"\n'
		antarray = AntArray(verbose=False)
		antarray.__dict__ = self.__dict__
		classhdf5 = jp.ClassHdf5(antarray, self.outname, verbose=False)
		classhdf5.Save()





	def Read( self, outname=None ) : 
		if (outname is not None) : self.outname = str(outname)
		if (self.outname[-5:] !='.hdf5') : self.outname += '.hdf5'
		if (self.verbose) : print 'AntArray.Read'
		if (self.verbose) : print '    Reading from  "'+jp.AbsPath(self.outname)+'"\n'
		antarray = AntArray(verbose=False)
		classhdf5 = jp.ClassHdf5(antarray, self.outname, verbose=False)
		classhdf5.Read()
		self.__dict__ = antarray.__dict__





	def Plotnfnv( self, nf, nv, returnlist=True ) : 
		'''
		nf: 
			index of antarray.freqorder
			nf = 0 ~ antarray.freqorder.size-1
			If nf==None, set nf=antarray.freqorder.size/2
		nv: 
			index of antarray.visorder, 
			nv = 0 ~ antarray.visorder.size-1
			If nf==None, set nv to the longest East-West baseline

		return:
			[nf, strfreq, nv, strblo, strbli, blname]
			nf: list of int
			strfreq: list of str like ['1400.0MHz', '1420.4MHz']
			nv: list of int
			strblo: list of str like ['2-3', '2-4']
			strblo: list of str like ['=(6.000, 0.000)', '=(10.400,5.987)'], has '=' in front of each element
			blname: one str
				cross: 'baseline='
				auto : 'channel='
				has '=' at the end

		returnlist:
			==True, nf, strfreq, nv, strblo, strbli are all list
			==False, return [nf[0], strfreq[0], nv[0], strblo[0], strbli[0], blname]
		'''
		# nf, strfreq are lists, antarray.vis[:,nf]
		try : nf = jp.npfmt(nf).flatten().round().astype(int)
		except : nf = None
		if (nf is None) : nf = jp.npfmt(self.vis.shape[1]/4)
		nf = nf[(0<=nf)*(nf<self.vis.shape[1])]
		strfreq = []
		for i in xrange(nf.size) : 
			try : strfreq.append(('%.1f' % self.freq[nf[i]])+'MHz')
			except : strfreq.append( 'axis1='+str(nf[i]) )
		#--------------------------------------------------

		# nv, strblo, strbli are lists, antarray.vis[:,:,nv]
		try : nv = jp.npfmt(nv).flatten().round().astype(int)
		except : nv = None
		vistype = 'cross' if(self.vis.dtype.name[:7]=='complex')else 'auto'
		blname = 'baseline=' if(vistype=='cross')else 'channel='
		strblo, strbli = [], []
		#--------------------------------------------------
		try : 
			bli = self.bli[:,:2]  # just use x,y
			if (nv is None) : nv = np.where(abs(bli[:,0])==abs(bli[:,0]).max())[0][:1]
			blo, bli = self.blo[nv], self.bli[nv]
			if (vistype == 'cross') : 
				for i in xrange(nv.size) : 
					if (self.Ant._telescope == 'paon4') : 
						n1, n2 = blo[i]/2 + 1
						if (blo[i,0]%2 == 1) : n1 = str(n1)+'H'
						else : n1 = str(n1)+'V'
						if (blo[i,1]%2 == 1) : n2 = str(n2)+'H'
						else : n2 = str(n2)+'V'
						strblo.append( n1+'-'+n2 )
					else : 
						strblo.append( str(blo[i,0])+'-'+str(blo[i,1])  )
					strbli.append( '=(%.3f, %.3f)' % tuple(bli[i]) )
			else : 
				for i in xrange(nv.size) : 
					if (self.Ant._telescope == 'paon4') : 
						n1, n2 = blo[i]/2 + 1
						if (blo[i,0]%2 == 1) : n1 = str(n1)+'H'
						else : n1 = str(n1)+'V'
						strblo.append( n1 )
					else : 
						strblo.append( str(blo[i,0]) )
					strbli.append( '' )
			#--------------------------------------------------
		except : 
			if (nv is None) : nv = jp.npfmt(self.vis.shape[2]/2)
			for i in xrange(nv.size) : 
				strblo.append( 'axis2='+str(nv[i]) )
				strbli.append( '' )
		#--------------------------------------------------

		if (not returnlist) : 
			nf, strfreq, nv, strblo, strbli = nf[0], strfreq[0], nv[0], strblo[0], strbli[0]
		return [nf, strfreq, nv, strblo, strbli, blname]



