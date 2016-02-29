#include <iostream>
#include <armadillo>
#include "main.h"
#include "dft_functions.h"

using namespace std;
using namespace arma;

extern "C" mat dualPotential(const operatorStruct Op, const paramStruct Pa,const mat dr,const cx_mat Sfi,const mat Xi,mat G2i)
{
  
  double omega=2.;
  arma::mat V(dr.n_rows,1);
  V=0.5*pow(omega,2)*(dr%dr);
  arma::mat Vdual=getVdual(Op,V);
  
  return Vdual;
  
}