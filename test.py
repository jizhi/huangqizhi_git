#! /usr/bin/env python
import numpy as np
import jizhipy as jp
import time


progressbar = jp.ProgressBar('Compledte:', 10)
for i in xrange(10) :
	progressbar.Progress()
	time.sleep(2)
	


