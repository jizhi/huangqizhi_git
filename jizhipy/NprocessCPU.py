
from ShellCmd import *


def NprocessCPU( Nprocess=None, warning=True ) : 
	uname = ShellCmd('uname')[0]
	if (uname == 'Linux') : 
		cores=ShellCmd('cat /proc/cpuinfo | grep "cpu cores"')[0]
		threads=ShellCmd('cat /proc/cpuinfo | grep siblings')[0]
	elif (uname == 'Darwin') : 
		cores = ShellCmd('sysctl machdep.cpu.core_count')[0]
		threads = ShellCmd('sysctl machdep.cpu.thread_count')[0]
	cores   = int(cores.split(':')[-1])
	threads = int(threads.split(':')[-1])
	cpuinfo = 'CPU INFO: '+str(cores)+' cores '+str(threads)+' threads'
	if (Nprocess is None): Nprocess = threads
	elif (Nprocess <= 0) : Nprocess = 1
	if (Nprocess>threads and warning) : Raise(Warning, cpuinfo+', but now  Nprocess='+str(Nprocess))
	return [Nprocess, cores, threads, cpuinfo]

