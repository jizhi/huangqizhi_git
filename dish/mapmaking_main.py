#! /usr/bin/env python
import jizhipy as jp
import numpy as np
from jizhipy.Plot import *

from AntArray import *
from Masking import *
from CaliGain import *
from CaliPhase import *
#from CaliTemperature import *
#from MapMaking import *
##################################################


print 'Start:', jp.Time(1)+'\n'
outdir = jp.Outdir(0, 'file')


#hdf5dir = '/disk/disk6/SunCygnusA_20160603190723_20160606080722/'
#pixstart, pixlength, pixperiod = 2324, 20, 2564-2324
#sourcename = 'CygA'
#brightsource = jp.BrightSource()


#hdf5dir = '~/study/lal/paon4/paon4_data/CygA/CygA665S1dec15_all-1_ac_0-4095-8'
hdf5dir = '/project/huangqizhi/paon4_data/CygA/CygA665S1dec15_all-1_ac_0-4095-8'
sourcename = 'CygA'
brightsource = jp.BrightSource()


##################################################


#
#antarraycross = AntArray(hdf5dir)
#
#antarraycross.WhichHdf5('transitsource')
#
#antarraycross.SelectVisType('cross1')
#antarraycross.SelectChannel([3, 7])
##antarray.SelectChannel([1, 3, 7])
#
##antarray.SelectVisType('auto1')
##antarray.SelectChannel([3])
#
##freq = 1400
##nfreq = abs(antarraycross.Ant.freq-freq)
##nfreq = np.where(nfreq==nfreq.min())[0][0]  # 307
##ntimemax = 69104
#
#
#
#maskingcross = Masking(antarraycross)
#
#axis, per, times = 0, 60, 1
##axis, per, times = 1, 7, 1
#nsigma, nloop, threshold = 4, 2, 0.001
#
##maskingcross.MaskLoop(axis=axis, per=per, times=times, nsigma=nsigma, nloop=nloop, threshold=threshold)
#
#
#
#caligaincross = CaliGain(antarraycross, maskingcross, outdir='cross1')
#
#caligaincross.Window(60*4)
#
#caligaincross.Gainnu()  # OK, very good
#
##caligain.See(3, 100, nfreq)
#
#masktime = (5000,11000)
#gainttimes = 100
#caligaincross.Gaint(masktime=masktime, gainttimes=gainttimes, legendsize=4, legendloc=8)
#
#
#
#nsigma = 4
#fwhm2sigmafactor = 2.287
#
##caliphase = CaliPhase(antarray, masking, caligain, nsigma=nsigma, fwhm2sigmafactor=fwhm2sigmafactor, Nprocess=200)
#caliphasecross = CaliPhase(antarraycross, None, None, nsigma=nsigma, fwhm2sigmafactor=fwhm2sigmafactor, Nprocess=200, plotfreq=1400, outdir='cross1')
#
#caliphasecross.RADec(sourcename)
#
#caliphasecross.Fringe()
#caliphasecross.Smooth(200)
#
#caliphasecross.FitBeam()
#
##caliphase.Ampb = np.load('mapmaking_main_output/CaliPhase_output/Ampb.npy')
##caliphase.Timeb = np.load('mapmaking_main_output/CaliPhase_output/Timeb.npy')
##caliphase.Sigmab = np.load('mapmaking_main_output/CaliPhase_output/Sigmab.npy')
##caliphase.Phasens = np.load('mapmaking_main_output/CaliPhase_output/Phasens.npy')
#
#caliphasecross.FitPhase()
#
##caliphase.Lewp = np.load('mapmaking_main_output/CaliPhase_output/Lewp.npy')
##caliphase.Phaseaddp = np.load('mapmaking_main_output/CaliPhase_output/Phaseaddp_4.npy')
#
##caliphasecross.Plot(dyDeff=0.5, dyLew=0.5, Nprocess=1)
##jp.Raise()
#


##################################################



antarrayauto = AntArray(hdf5dir)

antarrayauto.WhichHdf5('transitsource')

antarrayauto.SelectVisType('auto1')
antarrayauto.SelectChannel([3, 7])



maskingauto = Masking(antarrayauto)

axis, per, times = 0, 60, 1
#axis, per, times = 1, 7, 1
nsigma, nloop, threshold = 4, 2, 0.001

#maskingauto.MaskLoop(axis=axis, per=per, times=times, nsigma=nsigma, nloop=nloop, threshold=threshold)



caligainauto = CaliGain(antarrayauto, maskingauto, outdir='auto1')

caligainauto.Window(60*4)

caligainauto.Gainnu()  # OK, very good

#caligain.See(3, 100, nfreq)

masktime = (5000,11000)
gainttimes = 100
caligainauto.Gaint(masktime=masktime, gainttimes=gainttimes, legendsize=4, legendloc=8)



nsigma = 4
fwhm2sigmafactor = 2.287

#caliphase = CaliPhase(antarray, masking, caligain, nsigma=nsigma, fwhm2sigmafactor=fwhm2sigmafactor, Nprocess=200)
caliphaseauto = CaliPhase(antarrayauto, None, None, nsigma=nsigma, fwhm2sigmafactor=fwhm2sigmafactor, Nprocess=200, plotfreq=1400, outdir='auto1')

caliphaseauto.RADec(sourcename)

caliphaseauto.Fringe()
caliphaseauto.Smooth(200)
jp.Raise()



##################################################










calitemp = CaliTemp(antarray, antarray.Ant.freq.mean(), None)
calitemp.Vis()

pixstart, pixlength, pixperiod = masking.noisesource.pixstart, masking.noisesource.pixlength, masking.noisesource.pixperiod
calitemp.MaskNoiseSource(pixstart, pixlength, pixperiod)

maskmanual = np.zeros(calitemp.vis.shape, bool)
maskmanual[1500:] = True
calitemp.MaskManual(maskmanual)
maskmanual = 0 #@

calitemp.Gaint()
calitemp.vis.mask -= calitemp.masking.maskmanual
calitemp.vis.mask += calitemp.masking.masknoisesource

flux = brightsource.FluxDensity(sourcename, calitemp.freq)

calitemp.Smooth(201)
calitemp.CaliTemp(flux)
Raise()

calitemp.Smooth()
#calitemp.Plot()
calitemp.Smooth(201)

##################################################

mapmaking = MapMaking(caliphase, calitemp)
mapmaking.CaliPhase()
mapmaking.MapMaking(201)


print 'End:', jp.Time(1)
