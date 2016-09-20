import h5py
import sys
import numpy as np
import jizhipy as jp
##################################################


class _AntArrayBlorder( object ) : 


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


	def _Hdf5List( self, hdf5dir ) : 
		if (hdf5dir == '') : hdf5dir = './'
		if (hdf5dir[-1]!='/') : hdf5dir += '/'
		self.hdf5dir = hdf5dir
		self.hdf5list = jp.ShellCmd('ls '+hdf5dir+'*.hdf5')
		ndel = 0
		for i in xrange(len(self.hdf5list)) : 
			if (self.hdf5list[i-ndel][-9:] == '_old.hdf5') : 
				self.hdf5list.pop(i-ndel)
				ndel += 1


	def _WhichHdf5( self, whichhdf5 ) : 
		''' whichhdf5: can be absolute path of current hdf5(str), or the order of this hdf5(int)
		self.hdf5path: current hdf5 
		self.nhdf5: order of current hdf5 in hdf5list '''
		istype = jp.IsType()
		if (istype.isstr(whichhdf5)) : 
			self.hdf5path = whichhdf5
			self.nhdf5 = self.hdf5list.index(hdf5path)
		else : 
			self.nhdf5 = int(round(whichhdf5))
			self.hdf5path = self.hdf5list[self.nhdf5]


##################################################
##################################################
##################################################


class _AntArrayAnt( object ) : pass


##################################################
##################################################
##################################################


class AntArray( object ) : 
	dtype = 'class:'+sys._getframe().f_code.co_name
	''' This class just reads and handles the information in the HDF5 file, won't change anything ! 
	
	self.MaskChannel(), self.SelectChannel():
		Modify self.Blorder.maskorder, maskchannel, selectorder, selectchannel, auto, cross1, cross2

	self.SelectVisType():
		Modifies self.visorder

	self.Blorder.blorder, baseline, channelpos won't change anytime

	maskorder, selectorder, auto, cross1, cross2, visorder:
		Are used to self.vis[:,:,visorder], blorder[cross1], baseline[maskorder]

	maskchannel, selectchannel:
		Are used to channelpos[selectchannel]
	'''

	def __init__( self, hdf5dir=None ) : 
		'''
		self.__dict__.keys() 
			= [WhichHdf5(), MaskChannel(), SelectChannel(), SelectVisType(), Blorder, Hdf5, Ant, vis, vistype, visorder]

		self.Blorder.__dict__.keys() 
			blorder starts from 1, not 0 !
			= [blorder, auto1, auto2, cross1, cross2, cross3, feedpos, channelpos, baseline, Bl2Order(), Order2Bl(), maskchannel, maskorder, selectchannel, selectorder]
	
		self.Hdf5.__dict__.keys() 
			= [hdf5dir, hdf5list, hdf5path, nhdf5, obstime, sec1970,  nhdf5transit, transittimelocal, transittimetotal]

		self.Ant.__dict__.keys()
			= [noisesourcepos, inttime, lonlat, freq, dishdiam]

		hdf5dir:
			In this directory/folder, all .hdf5 are one continuum observation split into several files. If one file one observation, this hdf5dir must just contain one file.
		'''
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


	def WhichHdf5( self, whichhdf5 ) : 
		'''
		which:
			can be absolute path: which = '/xxx/yyy.hdf5'
			or index/number of this file: which = 3
			or 'transitsource' which is not the Sun
		'''
		istype = jp.IsType()
		if (istype.isstr(whichhdf5)) : 
			if (whichhdf5.lower() == 'transitsource') : which = 0
			else : which = whichhdf5
		else : which = whichhdf5
		self._WhichHdf5(which)
		if (istype.isstr(whichhdf5)) : 
			if (whichhdf5.lower() == 'transitsource') : 
				self._WhichHdf5(self.Hdf5.nhdf5transit[-1])


	def _WhichHdf5( self, whichhdf5 ) : 
		self.Hdf5._WhichHdf5(whichhdf5)
		fo = h5py.File(self.Hdf5.hdf5path, 'r')
		self.vis = fo['vis']  #@
		self.Hdf5.obstime = fo.attrs['obstime']
		self.Hdf5.sec1970 = fo.attrs['sec1970']
		# transittime (second) relates to the beginning of current hdf5 
		self.Hdf5.transittimelocal = fo['transitsource'][:].T[0] - self.Hdf5.sec1970
		# transittime (second) relates to the beginning of the first hdf5 (hdf5list[0])
		self.Hdf5.transittimetotal = self.Hdf5.transittimelocal + self.Hdf5.nhdf5*self.Ant.inttime*self.vis.shape[0]
		# nhdf5 of the transit source
		self.Hdf5.nhdf5transit = (self.Hdf5.transittimetotal / (self.Ant.inttime * self.vis.shape[0])).astype(int)


	def _Ant( self, fo ) : 
		self.Ant.lonlat = np.array([fo.attrs['sitelon'], fo.attrs['sitelat']])
		self.Ant.dishdiam = fo.attrs['dishdiam']
		self.Ant.inttime = fo.attrs['inttime']
		freq1 = np.arange(fo.attrs['freqstart'], fo.attrs['freqstart']+fo.attrs['freqstep']*fo.attrs['nfreq'], fo.attrs['freqstep'])
		freq2 = np.linspace(fo.attrs['freqstart'], fo.attrs['freqstart']+fo.attrs['freqstep']*fo.attrs['nfreq'], fo.attrs['nfreq'])
		self.Ant.freq = freq1 if(freq1.size==fo.attrs['nfreq'])else freq2
		self.Ant.noisesourcepos = np.array([-159.802, 11.920, 11.650])


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
		#--------------------------------------------------


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
			Raise(Warning, "vistype='"+vistype.lower()+"' not in ['auto1', 'auto2', 'cross1', 'cross2, cross3'], reset to vistype='cross1'")
			self.vistype = 'cross1'
			self.visorder = self.Blorder.cross1



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
