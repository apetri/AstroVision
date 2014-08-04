This python code makes use of the [LensTools](http://www.columbia.edu/~ap3020/LensTools/html/index.html) python package. To display the instructions just run

	python chi2.py -h

You can also enable multiprocessing via MPI

	mpiexec -n 3 chi2.py -f options.ini


Simulated noise maps
---------------------
You just need to run the script `generate_noise.py` in the following way

	python generate_noise.py 10
	
to generate for example 10 random noise realizations. After that you can load for example one of them in matlab with the following

	load(noise_map_1.mat)
	
	%The noise map is now available under the struct "noise map"
	
	size(noise_map)
	
	%should print [2048 2048]
	

MATLAB part
------------------------
getSimpleFeaure is for generating one-level features. Detailed description of the inputs are included in the function file itself.
		getSimpleFeature({'ref','om29','si850','w12','si750'},[-0.15,0.15],120,2.5,'scale25');

simpleChi is for computing the deltachi2.
