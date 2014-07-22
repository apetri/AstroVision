#Runs only in python3!
import sys,configparser

options = configparser.ConfigParser()
options["simulations"] = {"root_path":"/Users/andreapetri/Documents/Columbia/spurious_shear/convergence_maps"}
options["analysis"] = {
	
	"num_realizations" : "3",
	"smoothing_scales" : "0.1,0.5,1.0",
	"bin_edges" : "np.ogrid[-0.15:0.15:15j]",
	"redshift" : "1.0"

}

with open(sys.argv[1],"w") as configfile:
	options.write(configfile)