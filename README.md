## ksdft++
______

### In short:
	  ksdft++ simple code to compute electronic energies in atomic systems
	  by S.Konzett-Stoffl - 2015/16
	  after lecture by T.Arias (see http://dft.physics.cornell.edu/old-website/minicourse/dft.pdf)
	  or youtube for great lectures on dft {https://www.youtube.com/watch?v=oyvGeQ8ehBM}
	  
______
### External libraries used are:

	  for the matrix algebra:
	  
	  [armadillo](http://arma.sourceforge.net/) by Conrad Sanderson,
	  [released under MPL], see below:
	  (
	  Conrad Sanderson.
	  Armadillo: An Open Source C++ Linear Algebra Library for Fast Prototyping and Computationally Intensive Experiments.
	  Technical Report, NICTA, 2010. 
	  )
	  
	  for the fast fourier transform [fftw](http://www.fftw.org) 
	  
	  for input data processing ptree class of [boost](http://www.boost.org/) library is used
	 
______	 
### Features:
	  
	  as devised by T.Arias, formulation of matrix operations through
	  overloading of operators in operator classes
	  
	  thorough numerical testing routines, activate by setting
	  PERFORM_OPERATOR_TESTS in main.h
	  
	  automatic latex documentation, the document of which is generated 
	  in postprocessing routine plotScript.sh found in script folder

______
### Installation:

	  ksdft++ uses the make procedure; (aside you find the .cbp project file for processing with code::blocks ide)
	  ksdft++ only tested on ubuntu linux, may require adaptions to the parts using the system routines when porting to different platforms
	  
	  0.) pull the github repository:  
	      git pull https://github.com/sjjko/ksdft.git master
	      
	  
	  1.) install armadillo, boost, fftw libraries in their development versions:
	  sudo apt-get install libblas-dev
	  sudo apt-get install libboost-dev
	  sudo apt-get install libarmadillo-dev (need at least version 6.5)
	  
	  2.) adapt the library and include paths in Makfile to point to the directories where above libraries have been installed to
	  may look as follows:
	  LIB = -ldl /usr/local/lib/libarmadillo.so.6.500.5 /usr/lib/libfftw.so -lboost_iostreams -lboost_system -lboost_filesystem 
	  for the libraries to include and
	  INC = -I src/H -I /usr/local/include/armadillo_bits -I /usr/include 
	  for the include paths
	  
	  3.) do a bash setup.sh - should invoke make command and link script files
	  if make runs through without errors you should find binaries in /bin/Debug or/and /bin/Release respectively

______
### Running:

	  run ksdft++ by evoking script run.sh in script folder:
	  this script: 1.) cleans the directory
		       2.) runs ksdft++
		       3.) does postprocessing - generates images and converts latex output to pdf file

______
### Results:

	  After running you should find 
			  
			    - a pdf file documenting what ksdft++ has done in POSTPROCESSING/yourcasename/doc/yourcasename_doku.pdf
			    - images of wavefunctions and densities as computed in the POSTPROCESSING/yourcasename/IMG folder
			    - output of results data in POSTPROCESSING/yourcasename/DATA folder
			    
______
### Modification:

	  to add a new case proceed as follows:
	  
	  I) 	make new case directory in ksdft main directory case/caseName
	  II) 	copy a atoms.param file from Hatoms, which holds the atom coordinates, one in each row
	  III) 	copy a input.param file from Hatoms and change to your needs
	  IV) 	add an ionic potential in ionicPotentialClass::computePotential class method 
	  
	  

