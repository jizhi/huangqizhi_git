from Other import *
from npfmt import *
import multiprocessing
from ShellCmd import *
import signal





def NprocessCPU( Nprocess=None, verbose=True, hyper=True ) : 
	''' hyper: whether use hyper-threading, True/False '''
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
	if (Nprocess is None): 
		if (hyper) : Nprocess = threads
		else : Nprocess = cores
	elif (Nprocess <= 1) : Nprocess = 1
	if (Nprocess > threads and verbose) : Raise(Warning, cpuinfo+', but now  Nprocess='+str(Nprocess))
	return [Nprocess, cores, threads, cpuinfo]





class PoolFor( object ) : 
	'''
	def _DoMultiprocess( iterable ) : 
		return a
	
	pool = PoolFor(Nstart, Nend, Nprocess)
	data = pool.map_async(_DoMultiprocess, send, cast)
	
	data = np.concatenate(data, 0)
	'''


	def __init__( self, Nstart=None, Nend=None, Nprocess=None, nsplit=None, info=False, verbose=False ) :
		'''
		(1) PoolFor( Nstart, Nend, Nprocess )
				use Nstart, Nend, Nprocess to calculate nsplit
				split send in self.map_async()
		(2) PoolFor( nsplit )
				use this nsplit
				split send in self.map_async()
		(3) PoolFor()
				don't  split send in self.map_async(), send has been splitted when give to self.map_async( splitted_send, bcast )
				send[0] for process-1
				send[1] for process-2
				......
				send[n] for process-n+1
		'''
		self.zero, self.splitsend = False, True
		if ((Nstart is None or Nend is None) and nsplit is None) : 
			self.splitsend = False
			return
		#--------------------------------------------------
		if (nsplit is not None) : 
			Nprocess, cores, threads, cpuinfo = NprocessCPU(len(nsplit), verbose)
		#--------------------------------------------------
		else : 
			if (Nend-Nstart <= 0) : 
				Raise(Warning, 'Nend-Nstart='+str(Nend-Nstart)+'<=0, return None')
				self.zero = True
				return
			Nprocess, cores, threads, cpuinfo = NprocessCPU(Nprocess, verbose)
			if (Nend-Nstart < Nprocess) : Nprocess = 1
			# nsplit
			nsplit=np.linspace(Nstart,Nend, Nprocess+1).astype(int)
			nsplit = np.array([nsplit[:-1], nsplit[1:]]).T
		#--------------------------------------------------
		if (info) : print 'Open '+str(Nprocess)+' processes. '+cpuinfo
		self.Nprocess, self.nsplit = Nprocess, npfmt(nsplit)
		self.Nmax = (self.nsplit[:,1] - self.nsplit[:,0]).max()




	def map_async( self, func, send=None, bcast=None ) : 
		'''
		If use apply() or map(), we can stop the pool program with Ctrl-c. However, if use apply_async().get(xxx) or map_async().get(xxx), we can use Ctrl-c to stop the program at any time.

		iterable:
			iterable[0] = [n1, n2, PoolWorder-i], i=1,2,..,Nprocess

		func:
			_DoMultiprocess()

		send:
			None or tuple or 2D-ndarray
			If is tuple/list, means each element is one array (note that must be 2D, and split along axis=0/row)
			If not tuple/list, means the send is an unity: send = npfmt(send)
			If is 2D array, will split along axis-0 (row)

		bcast:
			None or tuple or others/as a unity
			If is others, may bcast = (bcast,)
			Some data which will be broadcasted to each processes, such as the x and p0 in test_multiprocess_poolfor_class-func.py
			Must be tuple: (x, p0, ...)
		'''
		if (self.zero) : return
		if (self.splitsend) : 
			istuple = True
			if (type(send) != tuple and type(send) != list) : 
				send, istuple = (npfmt(send),), False
			iterable, nsl = list(self.nsplit), self.nsplit
			for i in xrange(len(nsl)) : 
				iterable[i] = [ [tuple(iterable[i]), self.Nmax] ]
				sendtmp = ()
				for j in xrange(len(send)) : 
					if (send[j] is None) : sendtmp += (None,)
					else : 
						sendtmp += (send[j][nsl[i][0]:nsl[i][1]],)
				if (not istuple) : sendtmp = sendtmp[0]
				iterable[i] += [sendtmp, bcast]
		#--------------------------------------------------
		else : 
			self.Nprocess = len(send)
			iterable = []
			for i in xrange(len(send)) : 
				iterable.append( [[(None,None), self.Nmax], send[i], bcast] )
		#--------------------------------------------------
		pool = multiprocessing.Pool(self.Nprocess)
		self.data = pool.map_async(func, iterable).get(10**10)
		#--------------------------------------------------
		pool.close()
		pool.join()
		return self.data



