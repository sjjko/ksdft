#ifndef EXCVWN_H_INCLUDED
#define EXCVWN_H_INCLUDED

#include "main.h"
#include "customAssert.h"

using namespace arma;
using namespace myFunctions;

//arma::mat excVWN(arma::vec n);

// VWN parameterization of the exchange correlation energy
inline arma::mat excVWN(paramStruct Pa, arma::mat n)
{
// Constants:

const double X1 = 0.75*pow(3.0/(2.0*PI),2.0/3.0);
const double A = 0.0310907;
const double x0 = -0.10498;
const double b = 3.72744;
const double c = 12.9352;
const double Q = sqrt(4*c-b*b);
const double X0 = x0*x0+b*x0+c;
arma::mat rs(n),x(n),xones(n),X(n),out(n);
verbosity(Pa,"compute (4pi/3 n)^(⁻1/3)",2,__FILE__,__LINE__);
rs=pow(4*PI/3*n,-1/3); // Added internal conversion to rs
verbosity(Pa,"square root of result",2,__FILE__,__LINE__);
x=sqrt(rs);
xones=arma::ones<arma::mat>(x.n_rows,x.n_cols);
verbosity(Pa,"compute X vector values",2,__FILE__,__LINE__);
X=x%x+b*x+c*xones;
verbosity(Pa,"final assembly of exchange functional value into out matrix",2,__FILE__,__LINE__);
out=-X1/rs;
verbosity(Pa,"step 2",2,__FILE__,__LINE__);
out+= A*(log(x%x/X)+2*b/Q*atan(Q/(2*x+b*xones)));
verbosity(Pa,"step 3",2,__FILE__,__LINE__);
out-=A*((b*x0)/X0*(log((x-xones*x0)%(x-xones*x0)/X)+2*(2*x0+b)/Q*atan(Q/(2*x+b))));
verbosity(Pa,"sreturn result",2,__FILE__,__LINE__);
return out;
}

#endif // EXCVWN_H_INCLUDED
