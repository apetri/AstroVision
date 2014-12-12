from __future__ import print_function,division,with_statement
import argparse,ConfigParser

######################################################################
##################LensTools functionality#############################
######################################################################

from lenstools.simulations import IGS1
from lenstools.constraints import FisherAnalysis
from lenstools import Ensemble

###############################################################
####Borrow these from measure_all_features#####################
###############################################################

from measure_all_features import build_feature_list,Measurement

#############################################
#######Other functionality###################
#############################################

import astropy.units as u


################################################################################
#####This function computes the confusion matrix for a particular feature#######
################################################################################

def confusionMatrix(descriptor,measurement_list,measurement_covariance):

	#Instantiate a FisherAnalysis instance
	analysis = FisherAnalysis()

	#Populate with the models
	for measurement in measurement_list:
		ens = Ensemble.read(measurement.savename(descriptor))
		analysis.add_model(measurement.model.squeeze(with_ns=True),ens.mean())

	#Compute the covariance matrix
	covariance = Ensemble.read(measurement_covariance.savename(descriptor)).covariance()

	#Now we are ready to compute the confusion matrix, for each parameter
	return analysis,covariance




###############################################################
#########################Main##################################
###############################################################
def main(cmd_args):

	#Read options file
	options = ConfigParser.ConfigParser()
	options.read(cmd_args.options_file)
	save_path = options.get("analysis","save_path")
	smoothing_scale = options.getfloat("analysis","smoothing_scales") * u.arcmin

	#Build feature list
	feature_list = build_feature_list(options)

	#Construct all the IGS1 instances corresponding to the models we want to measure the features of
	all_igs1_models = IGS1.getModels(root_path=options.get("simulations","root_path"))

	#Build list for model book-keeping
	measurement_list = list()
	for model in all_igs1_models:
		measurement_list.append(Measurement(model=model,nrealizations=None,measurer=None,index=None,redshift=None,big_fiducial_set=False,smoothing=smoothing_scale,save_path=save_path))

	#Measurement of the covariance matrix
	measurement_covariance = Measurement(model=all_igs1_models[0],nrealizations=None,measurer=None,index=None,redshift=None,big_fiducial_set=True,smoothing=smoothing_scale,save_path=save_path)


	#Confusion matrix for first descriptor
	analysis,covariance = confusionMatrix(feature_list[0],measurement_list,measurement_covariance)

	return feature_list[0],analysis,covariance


if __name__=="__main__":

	#Parse command line options
	parser = argparse.ArgumentParser()
	parser.add_argument("-f","--file",dest="options_file",action="store",type=str,help="analysis options file")

	cmd_args = parser.parse_args()

	if cmd_args.options_file is None:
		parser.print_help()
		sys.exit(0)

	feature,analysis,covariance = main(cmd_args)

	
	


