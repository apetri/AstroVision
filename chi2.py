from __future__ import division,print_function,with_statement

import sys
import argparse,ConfigParser
import logging

####################################################
#######LensTools functionality######################
####################################################

from lenstools.simulations import IGS1

####################################################
###########Other functionality######################
####################################################

import numpy as np
from astropy.table import Table
from emcee.utils import MPIPool
from measure import measure_all_histograms

##########################################################
############Compute all the delta chi2####################
##########################################################

def compute_chi2(ensemble_list):
	return ensemble_list[0].compare(ensemble_list[1]),ensemble_list[0].compare(ensemble_list[2]),ensemble_list[0].compare(ensemble_list[3])

####################################################
#########Main#######################################
####################################################

if __name__=="__main__":

	#Initialize MPI pool
	try: 
		pool = MPIPool()
	except ValueError:
		pool = None

	#Parse command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-v","--verbose",action="store_true",default=False,dest="verbose",help="Display degug info")
	parser.add_argument("-f","--file",action="store",default="options.ini",dest="options_file",help="ini options file")
	parser.add_argument("-c","--compute",action="store_true",default=False,dest="compute",help="if enabled, computes the histograms, otherwise it just loads them from previously generated npy files")
	cmd_args = parser.parse_args()

	#Set logging level
	if cmd_args.verbose:
		logging.basicConfig(level=logging.DEBUG)
	else:
		logging.basicConfig(level=logging.INFO)

	if (pool is not None) and not(pool.is_master()):
	
		pool.wait()
		sys.exit(0)

	#Parse ini options file
	options = ConfigParser.ConfigParser()
	with open(cmd_args.options_file,"r") as configfile:
		options.readfp(configfile)

	#Build IGS1 instances for cosmological models handling
	fiducial_model = IGS1(name="fiducial",root_path=options.get("simulations","root_path"))
	high_Om_model = IGS1(name="high_Om",Om0=0.29,root_path=options.get("simulations","root_path"))
	low_w0_model = IGS1(name="low_w0",w0=-1.2,root_path=options.get("simulations","root_path")) 
	high_si8_model = IGS1(name="high_si8",sigma8=0.850,root_path=options.get("simulations","root_path"))

	models = [fiducial_model,high_Om_model,low_w0_model,high_si8_model]

	#Compute histogram ensembles for each of the models, otherways just load them from already npy generated file
	if cmd_args.compute:

		ensemble_array = measure_all_histograms(models,options,pool=pool)
		np.save("histograms.npy",ensemble_array)

	else:

		ensemble_array = np.load("histograms.npy")

	#Close pool if one is open
	if pool is not None:
		pool.close()

	#############################################################################################
	#######################Histograms are available here#########################################
	#######################in the ensemble_array structured array################################
	#############################################################################################

	#chi2 = compute_chi2(histogram_ensemble_list)
	
	#data_rows = [chi2]
	#t = Table(rows=data_rows,names=(r"$\Omega_m={0:.2f}$".format(high_Om_model.Om0),r"$w_0={0:.1f}$".format(low_w0_model.w0),r"$\sigma_8={0:.2f}$".format(high_si8_model.sigma8)))

	##################################
	#####Format the table#############
	##################################

	#for colname in t.columns:
	#	t[colname].format = "{0:.2f}"

	#t.write(sys.stdout,format="latex")

	logging.info("DONE!!")
