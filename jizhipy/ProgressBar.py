import sys
from Time import *
import numpy as np



class ProgressBar( object ) : 
	'''
	progressbar = ProgressBar()
	progressbar.__init__('RemoveRFI.RemoveLoop():', len(array))
	for i in xrange(len(a)) : progressbar.Progress()
	'''

	def __init__( self, string=None, Ntot=None, endnewline=True, pltend=False ):
		self.string, self.Ntot, self.count = string, Ntot, Ntot
		self.starttime0, self.starttime1 = Time()[:2]
		self.endnewline, self.pltend = endnewline, pltend
 

	def Progress( self, endstring='' ) : 
		self.count -=1
		currenttime0, currenttime1 = Time()[:2]
		dtime = Time(self.starttime0, currenttime0)
		if (self.count == self.Ntot-1) : font, end = '', ' \r'
		elif (self.count != 0) : font, end = '\r', ' '
		else : 
			if (self.endnewline or self.pltend) : font, end = '\r', ' \n'
			else : font, end = '\r', ' '
		string = font+self.string+'  '+str(self.Ntot-self.count)+'/'+str(self.Ntot)+'  '+dtime+endstring+end
		strblank = ' '*len(string)+'\r'
		sys.stdout.write(strblank)
		sys.stdout.flush()
		sys.stdout.write(string)
		sys.stdout.flush()
		if (self.pltend) : 
			if (self.count == 0) : 
				print self.starttime1+' => '+currenttime1

