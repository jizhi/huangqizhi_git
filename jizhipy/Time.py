import time
from IsType import *



def Time( a, b=None, c=None ) : 
	'''
	sec1970:
		second starts from '1970/01/01 00:00:00'
		It is an unity time all over the world (all timezones). No matter what timezone is, time.time() are the same !
	ephemtime:
		different "labels" of sec1970 at different timezones

	The time here is basing on the timezone, don't take local-law into account (Daylight Saving Time)!

	(1) Return local/system timezone (int)
	** Time('timezone')
	timezone >=0:east, <0:west

	(2) Convert ephemtime to sec1970 with timezone
	** Time(ephemtime, timezone)
	** Time('2016/10/25 16:39:24.68', 8)

	(3) Convert ephemtime from timezonein to timezoneout
	** Time(ephemtime, timezonein, timezoneout)
	** Time('2016/10/25 16:39:24.68', 8, 5)

	(4) Convert sec1970 to ephemtime with timezone
	** Time(sec1970, timezone)
	** Time(1477411200.729252, 8)

	(5) Calculate time interval between two ephemtimes
	** Time(ephemtime1, ephemtime2)
	** Time('2016/10/25 16:39:24.68', '2016/11/25 18:49:54.37')
	return (ephemtime2 - ephemtime1)

	(6) Return now-time with timezone WITH/WITHOUT Daylight Saving Time
	a in [0,1,2], b=timezone/None
	c=None/'daylight': use daylight, c=False: don't use daylight
	** Time(0/1/2, timezone, c)
	'''
	#--------------------------------------------------
	def Timezone( a, b ) : 
		'''
		Return local/system timezone (int)
		Time('timezone')
		>=0:east, <0:west
		'''
		return -time.timezone/3600
	#--------------------------------------------------
	def ephemTOsec1970( a, b ) : 
		'''
		ephemtime to sec1970 with timezone
		Time(ephemtime, timezone),  timezone can =None
		Time('2016/10/25 16:39:24.68', 8)
		'''
		if (b is None) : b, offset = Timezone('timezone',None), 0
		else : offset = b*3600 + time.timezone  # second
		# Calibrate the timezone without local-law (Daylight Saving Time)
		offset += time.mktime(time.strptime('1970/01/01 00:00:00', '%Y/%m/%d %H:%M:%S')) + Timezone('timezone',None)*3600
		n = a.rfind('.')
		if (n < 0) : dot = 0
		else : dot, a = float(a[n:]), a[:n]
		sec1970 = time.mktime(time.strptime(a, '%Y/%m/%d %H:%M:%S')) + dot - offset
		return sec1970
	#--------------------------------------------------
	def Interval( a, b ) : 
		'''
		time interval between two ephemtimes
		Time(ephemtime1, ephemtime2)
		ephemtime2 - ephemtime1
		Time('2016/10/25 16:39:24.68', '2016/11/25 18:49:54.37')
		'''
		n = a.rfind('.')
		if (n < 0) : dot = 0
		else : dot, a = float(a[n:]), a[:n]
		try : sec19701 = time.mktime(time.strptime(a, '%Y/%m/%d %H:%M:%S')) + dot 
		except : sec19701 = time.mktime(time.strptime(a, '%Y/%m/%d %p %I:%M:%S')) + dot 
		n = b.rfind('.')
		if (n < 0) : dot = 0
		else : dot, b = float(b[n:]), b[:n]
		try : sec19702 = time.mktime(time.strptime(b, '%Y/%m/%d %H:%M:%S')) + dot 
		except : sec19702 = time.mktime(time.strptime(b, '%Y/%m/%d %p %I:%M:%S')) + dot 
		d = (sec19702 - sec19701)/3600.  # hour
		h, d = int(d), (d-int(d))*60  # min
		m, s = int(d), (d-int(d))*60  # sec
		hms = str(h)+':'+str(m)+':'+('%.1f' % s)
		return hms
	#--------------------------------------------------
	def sec1970TOephem( a, b ) : 
		'''
		sec1970 to ephemtime with timezone
		Time(sec1970, timezone),  timezone can =None
		Time(1477411200.729252, 8)
		'''
		if (b is None) : b, offset = Timezone('timezone',None), 0
		else : offset = b*3600 + time.timezone  # second
		offset += time.mktime(time.strptime('1970/01/01 00:00:00', '%Y/%m/%d %H:%M:%S')) + Timezone('timezone',None)*3600
		dot, a = a-int(a), int(a)
		ephemtime = time.strftime('%Y/%m/%d %H:%M:%S', time.localtime(a+offset)) + str(dot)[1:4]
		return ephemtime
	#--------------------------------------------------
	def Nowtime( a, b, c ) : 
		'''
		Return now-time with timezone
		a in [0,1,2], b=timezone/None
		WITH 
		'''
		if (b is None) : b, offset = Timezone('timezone',None), 0
		else : offset = b*3600 + time.timezone  # second
		if (c is False) : offset += time.mktime(time.strptime('1970/01/01 00:00:00', '%Y/%m/%d %H:%M:%S')) + Timezone('timezone',None)*3600
		localtime = time.localtime(time.time()+offset)
		if (a == 0) : fmt = '%Y/%m/%d %H:%M:%S'
		elif (a == 1) : fmt = '%Y/%m/%d %p %I:%M:%S'
		elif (a == 2) : fmt = '%Y.%m.%d.%p.%I.%M.%S'
		timestr = time.strftime(fmt, localtime)
		return timestr
	#--------------------------------------------------
	istype = IsType()
	if (str(a).lower() == 'timezone') : 
		return Timezone(a, b)
	#--------------------------------------------------
	elif (istype.isstr(a) and (b is None or istype.isint(b))) : 
		if (c is None) : return ephemTOsec1970(a, b)
		else : return sec1970TOephem( ephemTOsec1970(a, b), c )
	#--------------------------------------------------
	elif (istype.isstr(a) and istype.isstr(b)) : 
		return Interval(a, b)
	#--------------------------------------------------
	elif (istype.isfloat(a) or a > 2) : 
		return sec1970TOephem(a, b)
	#--------------------------------------------------
	else : return Nowtime(a, b, c)


