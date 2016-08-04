


class ProgressBar( object ) : 
	'''
	progressbar = ProgressBar()
	progressbar.__init__('RemoveRFI.RemoveLoop():', len(array))
	for i in xrange(len(a)) : progressbar.Progress()
	'''

	def __init__( self, string=None, Ntot=None, pltend=False ):
		self.string, self.Ntot, self.count = string, Ntot, Ntot
		self.starttime0, self.starttime1 = Time()[:2]
		self.pltend = pltend
 
	def Progress( self ) : 
		self.count -=1
		currenttime0, currenttime1 = Time()[:2]
		dtime = Time(self.starttime0, currenttime0)
		if (self.count == self.Ntot-1) : font, end = '', '   \r'
		elif (self.count != 0) : font, end = '\r', '   '
		else : font, end = '\r', '   \n'
		string = font+self.string+'  '+str(self.Ntot-self.count)+'/'+str(self.Ntot)+'  '+dtime+end
		sys.stdout.write(string)
		sys.stdout.flush()
		if (self.pltend) : 
			if (self.count == 0) : 
				print self.starttime1+' => '+currenttime1

