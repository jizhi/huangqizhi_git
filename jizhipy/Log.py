import sys

'''
Import this module if you want to write a log file with argument "--log":
	python test.py --log
'''


if ('--log' in sys.argv[1:]) : 
	sysstd = SysStdout()
	sys.stdout = sysstd
	sys.stderr = sysstd



class LogStdout( object ) : 

	def __init__( self ) : 
		pyname = sys.argv[0].split('/')[-1]
		self.logf = open(pyname+'.log', 'a')
		self.logf.write('\n\n\n=============== '+Time(1)+' ===============\n')

	def write( self, outstream ) : 
		sys.__stdout__.write(outstream)
		if (outstream[-1]=='\r') : outstream = outstream[:-1]+'\n'
		if (outstream[0 ]=='\r') : 
			if (outstream[-1]=='\n') : outstream =outstream[1:]
			else : outstream = outstream[1:]+'\n'
		self.logf.write(outstream)
		self.logf.flush()

	def flush( self ) : 
		sys.__stdout__.flush()


