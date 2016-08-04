


def Time( time1=None, time2=None ):
	'''
	(1) time1==None and time2==None:
		return [timeephem, timesee, timeoutname]
	(2) time1 in [0,1,2] and time2==None
		return [timeephem, timesee, timeoutname][time1]
	(3) time1 is timeephem or timesee string, time2==None
		return sec1970
	(4) time1 and time2 are timeephem or timesee string
		return (time2 - time1) in second
	'''
	if (type(time1) != str) : 
		sec1970, localtime = time.time(), time.localtime()
		sec1970p = ('%.1f' % (sec1970-int(sec1970)))[1:]
		timeephem = time.strftime('%Y/%m/%d %H:%M:%S', localtime) + sec1970p
		timesee = time.strftime('%Y/%m/%d %p %I:%M:%S', localtime) + sec1970p
		timeoutname = time.strftime('%Y.%m.%d.%p.%I.%M.%S', localtime)
		timelist = [timeephem, timesee, timeoutname]
		if (time1 is None) : return timelist
		try : return timelist[time1]
		except : return timelist[0]
	else : 
		n = str(time1).rfind('.')
		if (n == -1) : n = len(time1)
		try : sec19701 = time.mktime(time.strptime(time1[:n], '%Y/%m/%d %p %I:%M:%S')) + float(time1[n:])
		except : sec19701 = time.mktime(time.strptime(time1[:n], '%Y/%m/%d %H:%M:%S')) + float(time1[n:])
		if (time2 is None) : return sec19701
		n = str(time2).rfind('.')
		if (n == -1) : n = len(time1)
		try : sec19702 = time.mktime(time.strptime(time2[:n], '%Y/%m/%d %p %I:%M:%S')) + float(time2[n:])
		except : sec19702 = time.mktime(time.strptime(time2[:n], '%Y/%m/%d %H:%M:%S')) + float(time2[n:])
		dt = (sec19702 - sec19701) /3600.
		h = int(dt)
		m = int((dt-h)*60)
		s = ((dt-h)*60-m)*60
		h, m, s = str(h), str(m), ('%.1f' % s)
		if (len(h) == 1) : h = '0'+h
		if (len(m) == 1) : m = '0'+m
		if (s.find('.') == 1) : s = '0'+s
		return h+':'+m+':'+s


