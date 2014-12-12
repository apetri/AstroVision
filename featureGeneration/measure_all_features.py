from __future__ import print_function,division,with_statement

import os,sys
import argparse,ConfigParser
import logging
import StringIO

######################################################################
##################LensTools functionality#############################
######################################################################

from lenstools import ConvergenceMap,GaussianNoiseGenerator
from lenstools.simulations import IGS1
from lenstools.index import Indexer,PowerSpectrum,PDF,Peaks,MinkowskiAll,Moments
from lenstools import Ensemble

########################################################################
##################Other functionality###################################
########################################################################

import numpy as np
from astropy.io import fits
from astropy.units import deg,arcmin
from emcee.utils import MPIPool

###########################################################################
#############Read INI options file and write summary information###########
###########################################################################

def write_info(options):

	s = StringIO.StringIO()

	s.write("""
Realizations to analyze: 1 to {0}

###########################################

""".format(options.get("analysis","num_realizations")))

	s.write("""Implemented descriptors
-------------------

""")

	if options.has_section("power_spectrum"):
		s.write("""Power spectrum: {0} bins between l={1} and l={2}\n\n""".format(options.get("power_spectrum","num_bins"),options.get("power_spectrum","lmin"),options.get("power_spectrum","lmax")))

	if options.has_section("moments"):
		s.write("""The set of 9 moments\n\n""")

	if options.has_section("peaks"):
		s.write("""Peak counts: {0} bins between kappa={1} and kappa={2}\n\n""".format(options.get("peaks","num_bins"),options.get("peaks","th_min"),options.get("peaks","th_max")))

	if options.has_section("minkowski_functionals"):
		s.write("""Minkowski functionals: {0} bins between kappa={1} and kappa={2}\n\n""".format(options.get("minkowski_functionals","num_bins"),options.get("minkowski_functionals","th_min"),options.get("minkowski_functionals","th_max")))

	s.seek(0)
	return s.read()

###########################################################################
########################IGS1 convergence maps measurer#####################
###########################################################################

def igs1_convergence_measure_all(realization,model,index,mask_filename=None,redshift=1.0,big_fiducial_set=False,smoothing=1.0*arcmin):

	"""
	Measures all the statistical descriptors of a convergence map as indicated by the index instance
	
	"""

	logging.debug("Processing {0}".format(model.getNames(realization,z=redshift,big_fiducial_set=big_fiducial_set,kind="convergence")))

	#Load the map
	conv_map = model.load(realization,z=redshift,big_fiducial_set=big_fiducial_set,kind="convergence")

	#Add the noise
	gen = GaussianNoiseGenerator.forMap(conv_map)
	noise = gen.getShapeNoise(z=redshift,ngal=15.0*arcmin**-2,seed=realization)

	logging.debug("Adding shape noise with rms {0:.3f}".format(noise.data.std()))
	conv_map += noise

	#Smooth the map
	logging.debug("Smoothing the map on {0}".format(smoothing))
	conv_map.smooth(scale_angle=smoothing)

	if mask_filename is not None:
		raise ValueError("Masks not implemented!") 
	
	logging.debug("Measuring...")

	#Allocate memory for observables
	descriptors = index
	observables = np.zeros(descriptors.size)

	#Measure descriptors as directed by input
	for n in range(descriptors.num_descriptors):

		
		if type(descriptors[n]) == PowerSpectrum:
			
			if mask_filename is None:
				l,observables[descriptors[n].first:descriptors[n].last] = conv_map.powerSpectrum(descriptors[n].l_edges)
			else:
				l,observables[descriptors[n].first:descriptors[n].last] = (conv_map*mask_profile).powerSpectrum(descriptors[n].l_edges)

		elif type(descriptors[n]) == Moments:

			if mask_filename is None:
				observables[descriptors[n].first:descriptors[n].last] = conv_map.moments(connected=descriptors[n].connected)
			else:
				observables[descriptors[n].first:descriptors[n].last] = masked_conv_map.moments(connected=descriptors[n].connected)
		
		elif type(descriptors[n]) == Peaks:
			
			if mask_filename is None:
				v,observables[descriptors[n].first:descriptors[n].last] = conv_map.peakCount(descriptors[n].thresholds,norm=descriptors[n].norm)
			else:
				v,observables[descriptors[n].first:descriptors[n].last] = masked_conv_map.peakCount(descriptors[n].thresholds,norm=descriptors[n].norm)

		elif type(descriptors[n]) == PDF:

			if mask_filename is None:
				v,observables[descriptors[n].first:descriptors[n].last] = conv_map.pdf(descriptors[n].thresholds,norm=descriptors[n].norm)
			else:
				v,observables[descriptors[n].first:descriptors[n].last] = masked_conv_map.pdf(descriptors[n].thresholds,norm=descriptors[n].norm)
		
		elif type(descriptors[n]) == MinkowskiAll:
			
			if mask_filename is None:
				v,V0,V1,V2 = conv_map.minkowskiFunctionals(descriptors[n].thresholds,norm=descriptors[n].norm)
			else:
				v,V0,V1,V2 = masked_conv_map.minkowskiFunctionals(descriptors[n].thresholds,norm=descriptors[n].norm)
			
			observables[descriptors[n].first:descriptors[n].last] = np.hstack((V0,V1,V2))
		
		elif type(descriptors[n]) == MinkowskiSingle:
			
			raise ValueError("Due to computational performance you have to measure all Minkowski functionals at once!")
		
		else:
			
			raise ValueError("Measurement of this descriptor not implemented!!!")

	#Return
	return observables


######################################################################################
##########Measurement object, handles the feature measurements from the maps##########
######################################################################################

class Measurement(object):

	"""
	Class handler for the maps feature measurements
	
	"""

	def __init__(self,model,nrealizations,measurer,index,redshift=1.0,big_fiducial_set=False,smoothing=1.0*arcmin,save_path=None):

		self.model = model
		self.nrealizations = nrealizations
		self.redshift = redshift
		self.measurer = measurer
		self.index = index
		self.big_fiducial_set = big_fiducial_set
		self.smoothing = smoothing

		#Build elements of save path for the features
		self.save_path = save_path

		#Cosmo id
		self.cosmo_id = self.model._cosmo_id_string
		if big_fiducial_set:
			self.cosmo_id += "_f"

		#Full save path
		dir_to_make = os.path.join(self.save_path,self.cosmo_id)
		if not os.path.isdir(dir_to_make):
			os.mkdir(dir_to_make)

		dir_to_make = os.path.join(dir_to_make,"smooth{0:02d}".format(int(self.smoothing.value*100)))
		if not os.path.isdir(dir_to_make):
			os.mkdir(dir_to_make)

		self.full_save_path = dir_to_make

	def savename(self,descriptor,save_type="npy"):

		return os.path.join(self.full_save_path,descriptor.name + ".{0}".format(save_type))		


	def measure(self,pool=None,save_type="npy"):
		"""
		Measures the features specified in the Indexer for all the maps whose names are calculated by get_all_map_names; saves the ensemble results in numpy array format

		"""

		realizations = range(1,self.nrealizations+1)

		#Build the ensemble
		ens = Ensemble.fromfilelist(realizations)

		#Load the data into the ensemble by calling the measurer on each map
		ens.load(callback_loader=self.measurer,pool=pool,model=self.model,index=self.index,mask_filename=None,redshift=self.redshift,big_fiducial_set=self.big_fiducial_set)

		#Break the ensemble into sub-ensemble, one for each feature
		single_feature_ensembles = ens.split(self.index)

		#For each of the sub_ensembles, save it in the appropriate directory
		for n,ensemble in enumerate(single_feature_ensembles):
			
			savename = self.savename(descriptor=self.index[n])
			logging.debug("Saving features to {0}".format(savename))
			ensemble.save(savename)


#######################################################
###############Main execution##########################
#######################################################

if __name__=="__main__":

	#Parse command line options
	parser = argparse.ArgumentParser()
	parser.add_argument("-f","--file",dest="options_file",action="store",type=str,help="analysis options file")
	parser.add_argument("-v","--verbose",dest="verbose",action="store_true",default=False,help="turn on verbosity")
	parser.add_argument("-t","--type",dest="type",action="store",default="npy",help="format in which to save the features")

	cmd_args = parser.parse_args()

	if cmd_args.options_file is None:
		parser.print_help()
		sys.exit(0)

	#Set verbosity level
	if cmd_args.verbose:
		logging.basicConfig(level=logging.DEBUG)
	else:
		logging.basicConfig(level=logging.INFO)

	#Initialize MPIPool
	try:
		pool = MPIPool()
	except:
		pool = None

	if (pool is not None) and not(pool.is_master()):
		
		pool.wait()
		sys.exit(0)
	
	logging.info("Start")

	#Parse INI options file
	options = ConfigParser.ConfigParser()
	with open(cmd_args.options_file,"r") as configfile:
		options.readfp(configfile)

	#Read the save path from options, along with the number of realizations, smoothing scale and redshift
	save_path = options.get("analysis","save_path")
	nrealizations = options.getint("analysis","num_realizations")
	smoothing_scale = options.getfloat("analysis","smoothing_scales") * arcmin
	redshift = options.getfloat("analysis","redshift")

	#Build an Indexer instance, that will contain info on all the features to measure, including binning, etc... (read from options)
	feature_list = list()

	if options.has_section("power_spectrum"):
		l_edges = np.ogrid[options.getfloat("power_spectrum","lmin"):options.getfloat("power_spectrum","lmax"):(options.getint("power_spectrum","num_bins")+1)*1j]
		np.savetxt(os.path.join(save_path,"ell.txt"),0.5*(l_edges[1:]+l_edges[:-1]))
		feature_list.append(PowerSpectrum(l_edges))

	if options.has_section("moments"):
		feature_list.append(Moments())

	if options.has_section("peaks"):
		th_peaks = np.ogrid[options.getfloat("peaks","th_min"):options.getfloat("peaks","th_max"):(options.getint("peaks","num_bins")+1)*1j]
		np.savetxt(os.path.join(save_path,"th_peaks.txt"),0.5*(th_peaks[1:]+th_peaks[:-1]))
		feature_list.append(Peaks(th_peaks))

	if options.has_section("minkowski_functionals"):
		th_minkowski = np.ogrid[options.getfloat("minkowski_functionals","th_min"):options.getfloat("minkowski_functionals","th_max"):(options.getint("minkowski_functionals","num_bins")+1)*1j]
		np.savetxt(os.path.join(save_path,"th_minkowski.txt"),0.5*(th_minkowski[1:]+th_minkowski[:-1]))
		feature_list.append(MinkowskiAll(th_minkowski))

	if options.has_section("pdf"):
		th_pdf = np.ogrid[options.getfloat("pdf","th_min"):options.getfloat("pdf","th_max"):(options.getint("pdf","num_bins")+1)*1j]
		np.savetxt(os.path.join(save_path,"th_pdf.txt"),0.5*(th_pdf[1:]+th_pdf[:-1]))
		feature_list.append(PDF(th_pdf))

	idx = Indexer.stack(feature_list)

	#Write an info file with all the analysis information
	with open(os.path.join(save_path,"INFO.txt"),"w") as infofile:
		infofile.write(write_info(options))


	#Construct all the IGS1 instances corresponding to the models we want to measure the features of
	all_igs1_models = IGS1.getModels(root_path=options.get("simulations","root_path"))

	fiducial_model = all_igs1_models[0]
	variation_models = all_igs1_models[1:]

	#First the fiducial model
	for big_fiducial_set in [True,False]:
		measurement = Measurement(model=fiducial_model,nrealizations=nrealizations,measurer=igs1_convergence_measure_all,index=idx,redshift=redshift,big_fiducial_set=big_fiducial_set,smoothing=smoothing_scale,save_path=save_path)
		logging.info("Processing {0}...".format(measurement.cosmo_id))
		measurement.measure(pool=pool,save_type=cmd_args.type)

	#Then all the others
	for model in variation_models:
		measurement = Measurement(model=model,nrealizations=nrealizations,measurer=igs1_convergence_measure_all,index=idx,redshift=redshift,big_fiducial_set=False,smoothing=smoothing_scale,save_path=save_path)
		logging.info("Processing {0}...".format(measurement.cosmo_id))
		measurement.measure(pool=pool,save_type=cmd_args.type)
	
	#Complete
	if pool is not None:
		pool.close()

	logging.info("DONE!")	

