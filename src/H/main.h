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

#undef PERFORM_OPERATOR_TESTS

#undef NORMALRUN
#undef FDTEST
#undef HATOMS

#ifdef HATOMS
    #undef SCHROEDINGER
    #undef FDTEST
#endif

//switches for energy components
    #undef CALC_KIN_ONLY
    #define CALC_KIN
    #define CALC_VEE
    #define CALC_EXC
    #define CALC_EION
    #define CALC_DEXC
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

#include "customAssert.h"

#define WITH_TEX
#undef USE_EXTERNAL_LIB_TEXCALLER // use external library texcaller - obsolete

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

#define PI arma::datum::pi

#define ISCRITICAL true
#define NOTCRITICAL false
//const bool ISCRITICAL=true;
//const bool NOTCRITICAL=true;

//arma::cx_vec func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward);
//arma::cx_vec func_dothefft(arma::cx_vec vec_Inp,arma::Col<double> S,double Faktor,fftw_direction dir);




#endif
