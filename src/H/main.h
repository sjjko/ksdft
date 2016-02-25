#ifndef MAIN_H
#define MAIN_H


// debug switches

#define CHECK_NAN_FINITE

// plot switches

#define PLOT_RESULT

#undef PLOT_PPM
#define PLOT_GNUPLOT

//switches for optimizers

#define DO_SD
#undef DO_PCCG

#define PERFORM_OPERATOR_TESTS

#undef NORMALRUN
#undef FDTEST
#define SCHROEDINGER
#undef HATOMS

#ifdef SCHROEDINGER
    #undef HATOMS
    #undef FDTEST
#endif
#ifdef HATOMS
    #undef SCHROEDINGER
    #undef FDTEST
#endif


//switches for energy components
#ifdef NORMALRUN
    #undef CALC_KIN_ONLY
    #define CALC_KIN
    #define CALC_VEE
    #define CALC_EXC
    #define CALC_EION
    #define CALC_DEXC
#endif
#ifdef SCHROEDINGER
    #undef CALC_KIN_ONLY
    #define CALC_KIN
    #undef CALC_VEE
    #undef CALC_EXC
    #define CALC_EION
    #undef CALC_DEXC
#endif
#ifdef FDTEST
    #undef CALC_KIN_ONLY
    #define CALC_KIN
    #define CALC_VEE
    #define CALC_EXC
    #define CALC_EION
    #define CALC_DEXC
#endif

#undef INCLUDE_PICTURE_IN_LATEX

#include <iostream>
#include <list>
#include <map>
#include <list>
#include <string>
#include <memory>
#include <exception>
//#include <fftw3.h>
#include <fftw.h>
#include <exception>
#include <armadillo>
#include "fft_minidft.h"
//#include "silo.h"
#include <iomanip>  // needed to use manipulators with parameters (precision, width)

extern int GLOBAL_verbosityLevel; //!< global variable set in main.cpp used in customAssert assert function
#include "customAssert.h"


#define WITH_TEX

#ifdef WITH_TEX
    #include <texcaller.h>
    //include to generate the image folder structure
    #include <sys/types.h>
    #include <sys/stat.h>
#endif

#define WITH_BOOST

#ifdef WITH_BOOST
    //!< includes if the BOOST parser mechanism is used
    #include <boost/foreach.hpp>
    #include <string>
    #include <set>
    #include <boost/property_tree/ptree.hpp>
    #include <boost/property_tree/ini_parser.hpp>
#endif

#define PLOT_WITH_CV //!< plot with cv_open library!

//#include "registryClass.h"
//#include "src/H/O.h"
//#include "L.h"
//#include "errorClass.h"

//#include "autoPtr.h"
//#include "error.h"
//#include "registryClass.h"
//#include "io.h"

#define PI arma::datum::pi

#define ISCRITICAL true
#define NOTCRITICAL false
//const bool ISCRITICAL=true;
//const bool NOTCRITICAL=true;

//arma::cx_vec func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward);
//arma::cx_vec func_dothefft(arma::cx_vec vec_Inp,arma::Col<double> S,double Faktor,fftw_direction dir);




#endif
