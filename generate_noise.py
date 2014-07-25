from __future__ import print_function

import sys

import numpy as np
import scipy.io as sio

if len(sys.argv)<2:
	print("Usage {0} <number_of_noise_realizations>".format(sys.argv[0]))
	sys.exit(1)

#Generate the white noise maps and save them
for n in range(1,int(sys.argv[1]) + 1):

	np.random.seed(n)
	noise_map = np.random.normal(loc=0.0,scale=0.47,size=(2048,2048))

	#Makes the noise map available through the "noise_map" struct in matlab
	sio.savemat("noise_map_{0}.mat".format(n),{"noise_map":noise_map})

#Done
print("Done!")
