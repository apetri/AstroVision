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

def compute_chi2_row(row):

	chi2_values = [row[1].compare(row[n]) for n in range(2,len(row))]
	return [row[0],] + chi2_values

def make_chi2_table(ensemble_array):
	
	#Compute all the delta chi2 and build table rows with them
	data_rows = [compute_chi2_row(ensemble_array[n]) for n in range(len(ensemble_array))]
	column_names = ("Smooth",) + ensemble_array.dtype.names[2:]

	#Build table
	return Table(rows=data_rows,names=column_names)


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

	#Close pool if one is open
	if pool is not None:
		pool.close()
	
	#Save computed histograms or load them
	if cmd_args.compute:
		np.save("histograms.npy",ensemble_array)
		logging.info("Saving histograms to histograms.npy")
	else:
		ensemble_array = np.load("histograms.npy")

	#############################################################################################
	#######################Histograms are available here#########################################
	#######################in the ensemble_array structured array################################
	#############################################################################################
	t = make_chi2_table(ensemble_array)

	t.write(sys.stdout,format="latex")

	logging.info("DONE!!")
