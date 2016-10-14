# Test suite for MOS-SIM

import numpy as np

def test():

	data = np.loadtxt('../data/W25000L1000', skiprows=2)

	print data

def rmsError():

	return