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

from emcee.utils import MPIPool
from measure import measure_all_histograms

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
	fiducial_model = IGS1(root_path=options.get("simulations","root_path"))
	high_Om_model = IGS1(Om0=0.29,root_path=options.get("simulations","root_path"))
	high_w0_model = IGS1(w0=-0.8,root_path=options.get("simulations","root_path")) 
	high_si8_model = IGS1(sigma8=0.850,root_path=options.get("simulations","root_path"))

	models = [fiducial_model]

	#Compute histogram ensembles for each of the models
	bin_edges,idx,histogram_ensemble_list = measure_all_histograms(models,options,pool=pool)

	#Close pool
	if pool is not None:
		pool.close()

	########################################################################################
	#######################Histograms are available here!!!#################################
	########################################################################################

	logging.info("DONE!!")