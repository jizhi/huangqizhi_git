from npfmt import *
from IsType import *


class BrightSource( object ) : 

	def __init__( self ) : 
		self.register_sources = [
		    'Cassiopeia_A','CasA', 
		    'Cygnus_A',    'CygA',
		    'Crab_Nebula', 'Crab', 'Taurus_A', 'TauA',
		    'Virgo_A',     'VirA']
	
		self.CasA=np.array([(23+23/60.+26/3600.)*15, 58+48/60.])
		self.CygA = np.array([(19+59/60.+28.4/3600)*15, +40+44/60.+2/3600.])
		self.Crab = self.TauA = np.array([(5+34/60.+31.9/3600)*15, +22+0+52.2/3600])
		self.VirA = np.array([(12+30/60.+49.4/3600)*15, +12+23/60+28/3600.])
	#	self.CygX = np.array([(20+31/60.)*15, +40+20/60.])


	def _CheckSourcename( self, sourcename ) : 
		if (sourcename not in self.register_sources) : Raise(Exception, "sourcename="+str(source)+" not in self.register_sources="+str(self.register_sources))


	def RADec( self, sourcename ) : 
		self._CheckSourcename(sourcename)
		sn = sourcename
		if   (sn in ['Cassiopeia_A','CasA']) : return self.CasA
		elif (sn in ['Cygnus_A','CygA']) : return self.CygA
		elif (sn in ['Crab_Nebula','Crab','Taurus_A','TauA']) : return self.Crab
		elif (sn in ['Virgo_A','VirA']) : return self.VirA


	def Convert( self, angle, what ) : 
		'''
		what: 'h2d' => hour to degree
		      'h2s' => hour to hms_str
		      'd2h' => degree to hour
            'd2s' => degree to hms_str
		      's2h' => hms_str to hour
            's2d' => hms_str to degree
		'''
		def h2s( ang ) : 
			h = int(ang[i])
			m = int((ang[i]-h)*60)
			s = ((ang[i]-h)*60-m)*60
			hms = str(h)+'h'+str(m)+'m'+('%.2f'%s)+'s'
			return hms
		def s2h( ang ) : 
			ang = str(ang)
			nh, nm = ang.find('h'), ang.find('m')
			h = float(ang[:nh])
			m = float(ang[nh+1:nm])
			s = float(ang[nm+1:-1])
			hms = h + m/60 + s/3600
			return hms
		#--------------------------------------------------
		istype = IsType()
		if (istype.isint(angle) or istype.isfloat(angle)) : 
			ang, islist = [angle], False
		else : ang, islist = angle, True
		ret = []
		for i in xrange(len(ang)) : 
			if   (what == 'h2d') : ret.append(ang[i]*15)
			elif (what == 'h2s') : ret.append(h2s(ang[i]))
			elif (what == 'd2h') : ret.append(ang[i]/15.)
			elif (what == 'd2s') : ret.append(h2s(ang[i]/15.))
			elif (what == 's2h') : ret.append(s2h(ang[i]))
			elif (what == 's2d') : ret.append(s2h(ang[i])*15)
		return ret


	def FluxDensity( self, sourcename, frequency ) : 
		'''
		sourcename: one str
	
		frequency: MHz, scale or ndarray
	
		Return:
			Flux density of this source at setting freq.
			Shape of return bases on shape of sourcename and freq.
			But each row of return is one source with several freq.

		Data from:
			The flux density values of standard sources used for antenna calibrations
				J.W.M. Baars, P.G. Mezger and H. Wendker
				1964
		'''
		freq = np.array([300.,400,408,500,600,700,750,800,900,
1000,1100,1200,1300,1400,1410,1420,1500,1600,1700,1800,1900,
2000,2500,2695,3000,3500,4000,4500,4995,5000,5500,6000,6500,
7000,7500,8000,8500,9000,9500,10000,10500,10690,11000,11500,
12000,12500,13000,13500,14000,14500,15000,15375,15500,16000,
16500,17000,17500,18000,18500,19000,19350,19500,20000,21000,
22000,23000,24000,25000,26000,27000,28000,29000,30000,31000,
31400,32000,33000,34000,35000,36000,37000,38000,39000,40000,
85000])
	
		flux_Cassiopeia_A = np.array([7609.,6145,6055,5206,
4546,4054,3852,3671,3364,3110,2898,2716,2560,2422,2410,2397,
2301,2194,2097,2010,1931,1858,1575,1489,1375,1226,1110,1017,
941,941,876,822,774,754,695,644,600,560,526,495,467,457,442,
420,399,380,363,347,333,319,307,298,295,284,274,265,256,247,
239,232,227,225,218,206,195,185,176,168,160,153,147,141,135,
130,128,125,121,117,113,109,106,102,99,96,40])
	
		flux_Cygnus_A = np.array([6375.,4948,4862,4065,3462,
3022,2844,2687,2422,2207,2029,1880,1752,1641,1631,1621,1544,
1459,1355,1266,1187,1116,856,783,689,574,489,425,376,375,335,
302,275,251,232,214,200,186,175,164,155,152,147,139,132,126,
120,115,110,106,102,99,98,94,91,87,84,82,79,77,75,74,72,68,
64,61,58,55,53,50,48,46,44,43,42,41,40,38,37,36,35,34,33,
32,13])
	
		flux_Crab_Nebula = np.array([1311,1221,1215,1156,1105,
1064,1046,1029,1000,974,951,931,913,896,895,893,881,867,854,
842,831,821,777,762,742,715,692,672,655,654,639,626,613,602,
592,583,574,566,558,551,545,542,539,533,527,522,517,512,507,
503,499,496,495,491,487,484,480,477,474,471,468,468,465,459,
454,449,444,440,436,431,428,424,420,417,416,414,411,408,405,
402,399,397,394,392,325])
	
		flux_Virgo_A = np.array([736.,583,573,486,419,369,349,
331,301,276,256,238,223,210,209,208,199,188,179,171,164,157,
131,123,113,100,89,81,75,75,69,64,60,57,54,51,48,46,44,42,
41,40,39,38,37,35,34,33,32,31,30,30,30,29,28,28,27,26,26,
25,25,25,24,23,22,22,21,20,19,19,18,18,17,17,17,16,16,16,15,
15,15,14,14,14,7])
	
		self._CheckSourcename(sourcename)
		sn = sourcename
		if   (sn in ['Cassiopeia_A','CasA']) : 
			flux = flux_Cassiopeia_A
		elif (sn in ['Cygnus_A','CygA']) : 
			flux = flux_Cygnus_A
		elif (sn in ['Crab_Nebula','Crab','Taurus_A','TauA']) : 
			flux = flux_Crab_Nebula
		elif (sn in ['Virgo_A','VirA']) : 
			flux = flux_Virgo_A
		istype = IsType()
		if (istype.isint(frequency) or istype.isfloat(frequency)) : islist = False
		else : islist = True
		frequency = npfmt(frequency, float)
		# power-law
		reflux = 10.**Interp1d( np.log10(freq), np.log10(flux), np.log10(frequency), 'linear' )
		if (not islist) : reflux = reflux[0]
		return reflux
	

