from Basic import *
import multiprocessing
from ShellCmd import *

##################################################



def NprocessCPU( Nprocess=None, warning=True ) : 
	try : Nprocess = int(round(Nprocess))
	except : Nprocess = None
	uname = ShellCmd('uname')[0]
	if (uname == 'Linux') : 
		cores=ShellCmd('cat /proc/cpuinfo | grep "cpu cores"')[0]
		threads=ShellCmd('cat /proc/cpuinfo | grep siblings')[0]
	elif (uname == 'Darwin') : 
		cores = ShellCmd('sysctl machdep.cpu.core_count')[0]
		threads = ShellCmd('sysctl machdep.cpu.thread_count')[0]
	cores   = int(cores.split(':')[-1])
	threads = int(threads.split(':')[-1])  # total
	cpuinfo = 'CPU INFO: '+str(cores)+' cores '+str(threads)+' threads'
	if (Nprocess is None): Nprocess = threads
	elif (Nprocess <= 1) : Nprocess = 1
	if (Nprocess>threads and warning) : Raise(Warning, cpuinfo+', but now  Nprocess='+str(Nprocess))
	return [Nprocess, cores, threads, cpuinfo]


##################################################
##################################################
##################################################


class PoolFor( object ) : 
	dtype = 'class:'+sys._getframe().f_code.co_name
	"""
	def _DoMultiprocess( iterable ) : 
		return a
	
	pool = PoolFor(Nstart, Nend, Nprocess)
	data = pool.map_async(_DoMultiprocess, send, cast)
	
	data = np.concatenate(data, 0)
	"""

	def __init__( self, Nstart, Nend, Nprocess=None, info=False, warning=True ) :
		self.zero = False
		if (Nend-Nstart <= 0) : 
			Raise(Warning, 'Nend-Nstart='+str(Nend-Nstart)+'<=0, return None')
			self.zero = True
			return
		Nprocess, cores, threads, cpuinfo = NprocessCPU(Nprocess, warning)
		if (Nend-Nstart < 2*Nprocess) : Nprocess = 1
		if (info) : print 'Open '+str(Nprocess)+' processes. '+cpuinfo
		self.Nstart, self.Nend, self.Nprocess = Nstart, Nend, Nprocess
		self.pool = multiprocessing.Pool(Nprocess)
		# nsplit
		nsplit = np.linspace(Nstart,Nend, Nprocess+1).astype(int)
		nsplit = list(np.array([nsplit[:-1], nsplit[1:]]).T)
		for i in xrange(len(nsplit)) : nsplit[i]=tuple(nsplit[i])
		self.nsplit = nsplit
		# self.nsplit = [(n1,n2), (n3,n4), (n5,n6), .....]


	def map_async( self, func, send=None, bcast=None ) : 
		'''
		If use apply() or map(), we can stop the pool program with Ctrl-c. However, if use apply_async().get(xxx) or map_async().get(xxx), we can use Ctrl-c to stop the program at any time.

		func:
			_DoMultiprocess()

		send:
			None or tuple or 2D-ndarray
			If is tuple, means each element is one array (note that must be 2D, and split along axis=0/row)
			If not tuple, means the send is an unity: send = npfmt(send)
			The whole 2D array which will be splited to each processes, such as the V in test_multiprocess_poolfor_class-func.py

		bcast:
			None or tuple or others/as a unity
			If is others, may bcast = (bcast,)
			Some data which will be broadcasted to each processes, such as the x and p0 in test_multiprocess_poolfor_class-func.py
			Must be tuple: (x, p0, ...)
		'''
		if (self.zero) : return
		if (type(bcast) == tuple) : 
			if (len(bcast) == 1) : bcast = bcast[0]
		if (type(send) != tuple) : send = (send,)
		iterable, nsplit = self.nsplit[:], self.nsplit[:]
		for i in xrange(len(nsplit)) : 
			iterable[i] = [iterable[i]]
			sendtmp = ()
			for j in xrange(len(send)) : 
				if (send[j] is None) : sendtmp += (None,)
				else : sendtmp += (send[j][nsplit[i][0]:nsplit[i][1]],)
			if (len(sendtmp) == 1) : sendtmp = sendtmp[0]
			iterable[i] += [sendtmp, bcast]
		#--------------------------------------------------
		self.data=self.pool.map_async(func, iterable).get(10**10)
		#--------------------------------------------------
		self.pool.close()
		self.pool.join()
		return self.data


