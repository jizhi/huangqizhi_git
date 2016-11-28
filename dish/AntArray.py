import h5py
import sys
import numpy as np
import jizhipy as jp
##################################################



class _AntArrayBlorder( object ) : 
	''' self.Blorder '''

	def _Blorder( self, blorder ) : 
		self.blorder = blorder[:]
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



class _AntArrayHdf5( object ) : 
	''' self.Hdf5 '''

	def _Hdf5List( self, hdf5dir ) : 
		''' self.hdf5dir, self.hdf5list '''
		if (hdf5dir == '') : hdf5dir = './'
		if (hdf5dir[-1] != '/') : hdf5dir += '/'
		self.hdf5dir = str(hdf5dir)
		self.hdf5list = jp.ShellCmd('ls '+hdf5dir+'*.hdf5')
		ndel = 0
		for i in xrange(len(self.hdf5list)) : 
			if (self.hdf5list[i-ndel][-9:] == '_old.hdf5') : 
				self.hdf5list.pop(i-ndel)
				ndel += 1



	def _WhichHdf5( self, which, Ant ) : 
		# self.Hdf5.nhdf5
		sumwhich = 1
		try : 
			which = np.sort(which).flatten()
			sumwhich = (which-np.arange(len(self.hdf5list))).sum()
		except : pass
		if (sumwhich == 0) : self.nhdf5 = which
		else : 
			for i in xrange(len(which)) : 
				try : which[i] = int(which[i])
				except : 
					try : which[i] = self.hdf5list.index(self.hdf5dir+which[i])
					except : which[i] = -1
			which = np.sort(np.array(which, int))
			self.nhdf5= which[(0<=which)*(which<len(self.hdf5list))]
			which= np.arange(self.nhdf5[0], self.nhdf5[-1]+1)
			if (self.nhdf5.size != which.size) : 
				missing = []
				for i in xrange(which.size) : 
					if (which[i] not in self.nhdf5) : missing.append(which[i])
				print 'Warning: self.Hdf5.nhdf5 is from '+str(which[0])+' to '+str(which[-1])+'(include), but now missing '+str(missing)+'\n'
		#--------------------------------------------------

		# self.Hdf5.hdf5path
		if (sumwhich != 0) : 
			self.hdf5path = []
			for i in xrange(self.nhdf5.size) : 
				self.hdf5path.append( self.hdf5list[self.nhdf5[i]] )
		else : self.hdf5path = self.hdf5list[:]
		#--------------------------------------------------

		fo = h5py.File(self.hdf5path[0], 'r')
		self.obstime = fo.attrs['obstime']
		self.sec1970 = fo.attrs['sec1970']
		#--------------------------------------------------

		# self.Hdf5.transittime (in second, start from hdf5list[0] file), self.Hdf5.nhdf5transit
		try : 
			transittime = fo['transitsource'][:].T[0]
			self.transittime = transittime - self.sec1970start
			self.nhdf5transit = (self.transittime / Ant.inttime / Ant.N0).astype(int)
		except : 
			if ('nhdf5transit' not in self.__dict__.keys()) : 
				print "Warning: NOT exist fo['transitsource'], you can set it by hand with AntArray.Transitsource()\n"
		#--------------------------------------------------



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
	''' This class just reads and handles the information in the HDF5 file, won't change anything ! 
	
	self.MaskChannel(), self.SelectChannel():
		Modify self.Blorder.maskorder, maskchannel, selectorder, selectchannel, auto1, auto2, cross1, cross2, cross3

	self.SelectVisType():
		Modifies self.visorder

	self.Blorder.blorder, baseline, channelpos won't change at any time !!!

	maskorder, selectorder, auto1, auto2, cross1, cross2, cross3, visorder:
		Are used to self.vis[:,:,visorder], blorder[cross1], baseline[maskorder]

	maskchannel, selectchannel:
		Are used to channelpos[selectchannel]
	'''


	def __init__( self, hdf5dir=None, verbose=True ) : 
		'''
		self.__dict__.keys() 
			= [WhichHdf5(), MaskChannel(), SelectChannel(), SelectVisType(), Blorder, Hdf5, Ant, vis, vistype, visorder]

		self.Blorder.__dict__.keys() 
			blorder starts from 1, not 0 !
			= [blorder, auto1, auto2, cross1, cross2, cross3, feedpos, channelpos, baseline, Bl2Order(), Order2Bl(), maskchannel, maskorder, selectchannel, selectorder]
	
		self.Hdf5.__dict__.keys() 
			= [hdf5dir, hdf5list, hdf5path, nhdf5, obstime, sec1970,  nhdf5transit, transittimelocal, transittimetotal]
				* transittimelocal, transittimetotal are pixels, NOT time in second

		self.Ant.__dict__.keys()
			= [noisesourcepos, inttime, lonlat, freq, antform]

		hdf5dir:
			In this directory/folder, all .hdf5 are one continuum observation splited into several files. If one file one observation, this hdf5dir must just contain one file.
		'''
		self.starttime = jp.Time(1)
		self.verbose = verbose
		if (self.verbose) : print '\n'
		if (hdf5dir is None) : return
		self.Blorder = _AntArrayBlorder()
		self.Hdf5    = _AntArrayHdf5()
		self.Ant     = _AntArrayAnt()
		self.Hdf5._Hdf5List(hdf5dir)
		fo = h5py.File(self.Hdf5.hdf5list[0], 'r')
		self._Ant(fo)
		self._Blorder(fo)
		self.MaskChannel()
		self.SelectVisType('all')
		self.Ant.N0 = fo['vis'].shape[0]
		self.Hdf5.sec1970start = fo.attrs['sec1970']
		self.Hdf5.obstimestart = fo.attrs['obstime']
		#--------------------------------------------------



	def WhichHdf5( self, whichhdf5 ) : 
		'''
		which:
			(1) Can be absolute path: which = '/xxx/yyy.hdf5'
			(2) Or index/number of the file: which = 3
				(1,2) whichhdf5 can be one or list: 
					whichhdf5 = [0,2,5,7,9,13,18]
					If that, will merge them into one array!
			(3) string with 'transitsource':
				1. =='transitsource': last transitsource
				2. =='transitsource5' or 'transitsource-5' or 'transitsource_5': 5-th transitsource (start from 1, not 0)
			(4) whichhdf5='all'
		'''
		istype = jp.IsType()
		if ('transitsource' in str(whichhdf5).lower()) : which = 0
		elif (str(whichhdf5).lower() == 'all') : which = np.arange(len(self.Hdf5.hdf5list))
		else : which = whichhdf5
		if (not (istype.islist(which) or istype.istuple(which) or istype.isndarray(which))) : which = [which]
		#--------------------------------------------------
		self.Hdf5._WhichHdf5(which, self.Ant)
		if ('transitsource' in str(whichhdf5).lower()) : 
				whichhdf5 = whichhdf5[len('transitsource'):]
				if (whichhdf5 == '') : nts = -1
				else : 
					whichhdf5 = whichhdf5.split('-')[-1]
					whichhdf5 = whichhdf5.split('_')[-1]
					nts = int(whichhdf5)
				try : nts = self.Hdf5.nhdf5transit[nts]
				except : jp.Raise(Exception, "NOT exist fo['transitsource'], you can set it by hand with AntArray.Transitsource()")
				self.Hdf5._WhichHdf5([nts], self.Ant)
		#--------------------------------------------------



	def Transitsource( self, transittime ) : 
		''' 
		If NO fo['transitsource'], set self.Hdf5.nhdf5transit and self.Hdf5.transittime (second from hdf5list[0]) here

		transittime:
			in second, start from hdf5list[0] (the beginning)
		'''
		self.Hdf5.transittime = jp.npfmt(transittime).astype(float).take(0)
		self.Hdf5.nhdf5transit = self.Hdf5.transittime / self.Ant.inttime / self.Ant.N0
		#--------------------------------------------------



	def _Ant( self, fo ) : 
		self.Ant.lonlat = np.array([fo.attrs['sitelon'], fo.attrs['sitelat']])
		self.Ant.inttime = fo.attrs['inttime']
		self.Ant.telescope = fo.attrs['telescope']
		try : self.Ant.noisesourcepos = fo['noisesource'][0]
		except : pass
		#--------------------------------------------------
		if ('cylinder' in self.Ant.telescope.lower()) : 
			self.Ant.antform = np.array([fo.attrs['cywid'], fo.attrs['cylen']])
		else : self.Ant.antform = fo.attrs['dishdiam']
		#--------------------------------------------------
		freq1 = np.arange(fo.attrs['freqstart'], fo.attrs['freqstart']+fo.attrs['freqstep']*fo.attrs['nfreq'], fo.attrs['freqstep'])
		freq2 = np.linspace(fo.attrs['freqstart'], fo.attrs['freqstart']+fo.attrs['freqstep']*fo.attrs['nfreq'], fo.attrs['nfreq'])
		self.Ant.freq = freq1 if(freq1.size==fo.attrs['nfreq'])else freq2
		#--------------------------------------------------
		self.freqorder = np.arange(self.Ant.freq.size)



	def _Blorder( self, fo ) : 
		self.Blorder._Blorder(fo['blorder'])
		self.Blorder._AutoCross()
		self.Blorder._Feedpos(fo['feedpos'])
		self.Blorder._Baseline()



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
		_selectorder=np.ones(len(self.Blorder._blorder),bool)
		_selectorder[self.Blorder.maskorder] = False
		self.Blorder.selectorder = np.arange(_selectorder.size)[_selectorder]
		#--------------------------------------------------
		self._MaskChannel()  # Recal auto, cross1, cross2
		self.SelectVisType(self.vistype)  # Reset visorder



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
			if (self.verbose) : print "Warning: vistype='"+vistype.lower()+"' not in ['auto1', 'auto2', 'cross1', 'cross2, cross3'], reset to vistype='cross1'\n"
			self.vistype = 'cross1'
			self.visorder = self.Blorder.cross1
		if (self.visorder.size == 0) : jp.Raise(Exception, 'self.visorder.size == 0')



	def SelectFreq( self, nfreq, unit=None ) : 
		'''
		nfreq, unit:
			If ('hz' not in str(unit).lower()) : 
					self.freqorder = nfreq
					antarray.vis[:, self.freqorder, self.visorder]
			Else, use nfreq to get self.freqorder
		'''
		nfreq = jp.npfmt(nfreq)
		if ('hz' not in str(unit).lower()) : unit = ''
		else : unit = unit.lower()
		if (unit == '') : self.freqorder = nfreq.astype(int)
		else : 
			if (unit == 'hz') : nfreq = nfreq * 1e-6  # MHz
			elif (unit == 'ghz') : nfreq = nfreq * 1e3
			else : nfreq = 1.*nfreq
			dfreq = abs(nfreq[:,None] - self.Ant.freq[None,:])
			self.freqorder= np.where(dfreq==dfreq.min(1)[:,None])[1]
		#--------------------------------------------------



	def Vis( self ) : 
		'''
		Read the visibibity data

		self.vis:
			np.concatenate( fo['vis'][:,freqorder,visorder] )

		self.timem:
			in minute
			self.timem[0] is the time(minute) on day self.Hdf5.obstime.split(' ')[0] (start from 0:0:0 A.M.)
			Note that it is local time!
		'''
		if (self.verbose) : 
			starttime = jp.Time(1)
			print 'AntArray.Vis: start @', starttime

		self.vis, axis = [], 1
		if (self.freqorder.size > self.visorder.size) : axis = 2

		for i in xrange(len(self.Hdf5.hdf5path)) : 
			if (self.verbose) : print '    Reading file '+str(self.Hdf5.nhdf5[i])+' = "'+self.Hdf5.hdf5path[i].split('/')[-1]+'"'
			fo = h5py.File(self.Hdf5.hdf5path[i], 'r')
			vistmp = []
			if (axis == 1) : 
				for j in xrange(self.freqorder.size) : 
					if ('auto' in self.vistype) : vistmp.append( fo['vis'][:,self.freqorder[j],self.visorder][:,None].real )
					else : vistmp.append( fo['vis'][:,self.freqorder[j],self.visorder][:,None] )
				vistmp = np.concatenate(vistmp, 1)
				if (len(vistmp.shape)==2) : vistmp = vistmp[:,:,None]
			elif (axis == 2) : 
				for j in xrange(self.visorder.size) : 
					if ('auto' in self.vistype) : vistmp.append( fo['vis'][:,self.freqorder,self.visorder[j]][:,:,None].real )
					else : vistmp.append( fo['vis'][:,self.freqorder,self.visorder[j]][:,:,None] )
				vistmp = np.concatenate(vistmp, 2)
			self.vis.append(vistmp)
		self.vis = np.concatenate(self.vis, 0)
		#--------------------------------------------------

		localtime = np.array(self.Hdf5.obstime.split(' ')[-1].split(':'), float)
		localtime = localtime[0]*60+localtime[1]+localtime[2]/60
		self.timem = np.arange(len(self.vis)) * self.Ant.inttime /60. + localtime 
		#--------------------------------------------------
		if (self.verbose) : 
			endtime = jp.Time(1)
			costtime = jp.Time(starttime, endtime)
			tottime = jp.Time(self.starttime, endtime)
			print 'AntArray.Vis:  end  @', endtime
			print 'Cost:', costtime+'     Total:', tottime+'\n'

		










	#---------- def _Ant() ----------
	def LonLat( self, lon, lat ) : 
		''' longitude and latitude of the antenna, in degree '''
		self.Ant.lonlat = np.array([lon, lat]).flatten()

	def Diameter( self, diameter ) : 
		''' Diameter of the single antenna, in meter '''
		self.Ant.dishdiam = np.array(diameter).take(0)

	def Inttime( self, inttime ) : 
		''' Integration time, in second '''
		self.Ant.inttime = np.array(inttime).take(0)

	def Freq( self, freq ) : 
		''' freqlist, in MHz. 1D '''
		self.Ant.freq = jp.npfmt(freq).flatten()

	def Noisesourcepos( self, x, y, z ) : 
		''' Coordinate of the noise source on the hill, meter '''
		self.Ant.noisesourcepos = np.array([x, y, z]).flatten()
	#---------- def _Ant()  END ----------


	#---------- def _Blorder() ----------
	def BlorderSet( self, blorder ) : 
		self.Blorder._Blorder(blorder)
		self.Blorder._AutoCross()
		keys = self.Blorder.__dict__.keys()
		if ('blorder' in keys and 'feedpos' in keys) : 
			self.MaskChannel()
			self.SelectVisType('all')

	def FeedposSet( self, feedpos ) : 
		''' Feeds position, in meter
		Order: feed-1,  feed-2,  feed-3,  feed-4, ...,  feed-n
		       Ant-1-H, Ant-1-V, Ant-2-H, Ant-2-V, ..., Ant-n/2-V
		'''
		self.Blorder._Feedpos(feedpos)
		self.Blorder._Baseline()
		keys = self.Blorder.__dict__.keys()
		if ('blorder' in keys and 'feedpos' in keys) : 
			self.MaskChannel()
			self.SelectVisType('all')
	#---------- def _Blorder()  END ----------
