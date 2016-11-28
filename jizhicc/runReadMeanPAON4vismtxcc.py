#! /usr/bin/env python
import os, sys
##################################################
##################################################
##################################################




# ReadMeanPAON4vismtx.cc is the c++ program to read, merge and average the .ppf files
# Qizhi Huang, 23 November 2015



# If we modified ReadMeanPAON4vismtx.cc, you need to compile it before run it. Default is False
#recompile = True   # False
recompile = False



# PAON4 has 36 visibilities. Valid visibilities: 8 auto-correlations, 6 cross-correlations * 2 polarizations, total 8+6*2=20.
# acORall = 'all' : save all 36 visibilities.
# acORall = 'ac' : just save 20 vaild visibilities to decrease the size of output file.
#acORall = 'all'
acORall = 'ac'



# How many .ppf files you want to read and merge?
# times = 'x,y,z,n' : read and merge .ppf files from 'x'(include) to 'y'(include), then average them every 'z' files, then save it into n parts.
# If you want to read and merge all .ppf files of one data set (for example, CygA665S1dec15), but you don't know how many files in it, you can set y='max'.
# Example (1) time = '0,max,1' or '0,max,1,1' : read and merge all .ppf files, don't average them (z=1), save it to one file.
# Example (2) time = '100,5000,400,3' : read and merge .ppf files from 100(include) to 5000(include), average them every 400 files, then save it into 3 filts/parts.
# If the data with time exceeds the memory of the machine, then devides it into smaller pieces automatically. So, it's convenient to set time='0,max,1'
time = '0,max,1'



# Similar as time above, but without ",n".
# Example (1) freq = '0,4095,8' : save frequency bins from 0(include) ot 4095(include), and average it every 8 frequency bins.
# Example (2) you can also just save frequency around 21cm (1416-1425MHz) : freq = '2730,2868,1'.
# Average over the frequency can decrease the size of output file.
freq = '0,4095,1'
#freq = '0,4095,8'



# Which folder to save the output files.
# outdir = '' or './' : save to current folder.
# outdir = '../' : save to the last level folder.
# outdir = '../haha/' : save to folder '../haha/'
outdir = 'MergePPF/'



# Which data set you want to read and merge.
# (1) can be string : ppfdirlist = 'CygA665S1dec15'
# (2) can be list of string : ppfdirlist = ['CygA565S2dec15', 'CygA665S1dec15', 'CygA765S30nov15']
#ppfdirlist = ['CygA6sep16', 'CasA26sep16', 'CygA14sep16', 'CasA29sep16', 'CasA30sep16']
ppfdirlist = ['CygA6sep16', 'CygA14sep16'][:1]
#ppfdirlist = 'CasA8mai16'



##################################################
##################################################
##################################################



if (outdir != '') : 
	if (outdir[-1] != '/') : outdir += '/'
if (outdir not in ['', './', '../']) : 
	outdir = os.path.expanduser(outdir)
	if (not os.path.exists(outdir)) : os.mkdir(outdir)


times = time.split(',')
if (len(times) == 3) : Npart = 1
else : Npart = int(times[3])
t0, dt = int(times[0]), int(times[2])
try : t1 = int(times[1])
except : t1 = 'max'


freqs = freq.split(',')
strf = freqs[0]+'-'+freqs[1]+'-'+freqs[2]
f0, f1, df = int(freqs[0]), int(freqs[1]), int(freqs[2])


codename = 'ReadMeanPAON4vismtx.cc'
if (type(ppfdirlist) == str) : ppfdirlist = [ppfdirlist]
if (ppfdirlist == []) : 
	print 'ppfdirlist = []'
	sys.exit()


#--------------------------------------------------
''' For my mac '''
#command1 = 'c++ -g -DDarwin -I/usr/local/Sophya/include/ -I/opt/local/include -I/usr/X11R6/include/ -fno-common -O -fPIC -c -I./inlib -c -o '+codedir+'a.o '+codedir+codename
#command2 = 'c++ -fno-common -g -O -fPIC -bind_at_load -L/usr/local/Sophya/slb/ -lPI -lextsophya -lsophya -L/opt/local/lib -lXm -ledit -lcurses -L/usr/X11R6/lib/ -lXt -lX11 -L/usr/local/Sophya/ExtLibs/lib -lcfitsio -lfftw3 -lfftw3f -lfftw3_threads -lfftw3f_threads -lxastro -framework Accelerate -lpthread -lm -lc -ldl -o '+codedir+codename[:-3]+'.out '+codedir+'a.o'


''' For bao@bao3 '''
CXXCOMPLIE = 'g++ -DLinux -I/ope/local/include -I/Dev/Sophya64/include -I/Dev/ExtLibs/include -I/usr/X11R6/include -Wall -Wpointer-arith -fno-common -O -g -fPIC -c '
CXXLINK = 'g++ -Wall -Wpointer-arith -O -g -fPIC '
SOPHYAEXTSLBLIST = '-L/Dev/Sophya64/slb -lextsophya -lsophya -lPI -L/Dev/ExtLibs/lib -lcfitsio -lfftw3 -lfftw3f -lfftw3_threads -lfftw3f_threads -llapack -lblas -lxastro -lgfortran -lstdc++ -lpthread -lm -lc -ldl '
command1 = CXXCOMPLIE+'-o a.o '+codename
command2 = CXXLINK+SOPHYAEXTSLBLIST+'-o '+codename[:-3]+'.out a.o '
#--------------------------------------------------


if (recompile) : 
	if (os.path.exists(codename[:-3]+'.out')) : 
		os.system('rm '+codename[:-3]+'.out')
	os.system(command1)
	os.system(command2)
	os.system('rm a.o')



##################################################
##################################################
##################################################



sysmem = os.popen('free -m').readlines()[1]
for i in xrange(len(sysmem)) : 
	if (sysmem[i] in '0123456789') : break
for j in xrange(i, len(sysmem)) : 
	if (sysmem[j] == ' ') : break
sysmem = float(sysmem[i:j])*1



for i in range(len(ppfdirlist)) : 
	ppfdir = ppfdirlist[i]


	if (t1 == 'max') : # include t1
		t1 = (os.popen('ls /Raid-bao5/PAON4/'+ppfdir+' | wc').readlines())[0].split(' ')
		for j in xrange(len(t1)) : 
			if (t1[j] != '') : 
				t1 = int(t1[j])/4-2
				break

	# time=t0,t1,dt,Npart
	# freq=f0,f1,df
	num = 1.5e8
	if (acORall == 'all') : nv = 36 *2  # (*2 because of complex)
	elif (acORall == 'ac') : nv = 2 * (6*2 + 4)
	arrmem = 1.*(t1-t0)/dt/Npart *(f1-f0)/df *nv /num *800 # MB

	if (arrmem > sysmem) : 
		arrmem = 1.*(t1-t0)/dt *(f1-f0)/df *nv /num *800
		Npart = arrmem / sysmem
		if (Npart == int(Npart)) : Npart = int(Npart)
		else : Npart = int(Npart)+1


	Ntime = (1.+t1) /Npart /dt
	if (Ntime == int(Ntime)) : Ntime = int(Ntime)
	else : Ntime = int(Ntime)+1

	tlist = [0]
	for j in xrange(Npart) : 
		tlist.append( tlist[-1] + dt*Ntime )
	for j in xrange(len(tlist)) : 
		if (tlist[j] > t1+1) : k = j+1
	try : tlist = tlist[:k]
	except : pass
	if (tlist[-2] == t1+1) : tlist = tlist[:-1]
	else : tlist[-1] = t1+1

	strt, time = [], []
	for j in xrange(len(tlist)-1) : 
		strt.append(str(tlist[j])+'-'+str(tlist[j+1]-1)+'-'+str(dt))
		time.append(str(tlist[j])+','+str(tlist[j+1]-1)+','+str(dt))


	for j in xrange(len(strt)) : 
		outname = outdir+ppfdir+'_'+strt[j]+'_'+acORall+'_'+strf+'.fits'
		cmd = './ReadMeanPAON4vismtx.out /Raid-bao5/PAON4/'+ppfdir+' /Raid-bao6/PAON4/'+ppfdir+' '+time[j]+' '+acORall+' '+freq+' '+outname
		print cmd
		os.system(cmd)
		sys.exit()


