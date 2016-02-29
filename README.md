ksdft++ simple code to compute electronic energies in atomic systems
by S.Konzett-Stoffl - 2015/16
after lecture by T.Arias (see http://dft.physics.cornell.edu/old-website/minicourse/dft.pdf)
or youtube for great lectures on dft {https://www.youtube.com/watch?v=oyvGeQ8ehBM}

External libraries used are:

	  for the matrix algebra:
	  
	  armadillo (http://arma.sourceforge.net/) by Conrad Sanderson,
	  [released under MPL], see below:
	  (
	  Conrad Sanderson.
	  Armadillo: An Open Source C++ Linear Algebra Library for Fast Prototyping and Computationally Intensive Experiments.
	  Technical Report, NICTA, 2010. 
	  )
	  
	  for the fast fourier transform:
	  
	  fftw (http://www.fftw.org) 
	  [released under the GNU General Public License]
	  
	  for input data processing:
	  
	  libboost: ptree
	  
Features:
	  
	  as devised by T.Arias, formulation of matrix operations through
	  overloading of operators in operator classes
	  
	  thorough numerical testing routines, activate by setting
	  PERFORM_OPERATOR_TESTS in main.h
	  
	  automatic latex documentation, the document of which is generated 
	  in postprocessing routine plotScript.sh found in script folder
	  
Installation:

	  ksdft++ uses the make procedure; (aside you find the .cbp project file for processing with code::blocks ide)
	  
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
	  
	  3.) make debug OR make release
	  if make runs through without errors you should find binaries in /bin/Debug or/and /bin/Release respectively
	  
	  4.) make a soft link to script run.sh in script directory
	  
	  5.) run ksdft++ by evoking script run.sh 
	  this script: 1.) cleans the directory
		       2.) runs ksdft++
		       3.) does postprocessing - generates images and converts latex output to pdf file
	  
Modification:

	  to add a new case proceed as follows:
	  
	  make new case directory in ksdft main directory case/caseName
	  copy a atoms.param file from Hatoms, which holds the atom coordinates, one in each row
	  copy a input.param file from Hatoms and change to your needs
	  add an ionic potential in ionicPotentialClass::computePotential class method 
	  
	  

