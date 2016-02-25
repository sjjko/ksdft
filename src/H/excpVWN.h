#ifndef EXCPVWN_H_INCLUDED
#define EXCPVWN_H_INCLUDED

#include "main.h"

#include <armadillo>
#include "customAssert.h"


// VWN parameterization of the exchange correlation energy
inline arma::mat excpVWN(arma::mat n)
{

//function out=excpVWN(n)
// Constants
const double X1 = 0.75*pow((3.0/(2.0*PI)),(2.0/3.0));
const double A = 0.0310907;
const double x0 = -0.10498;
const double b = 3.72744;
const double c = 12.9352;
double Q = sqrt(4.*c-b*b);
double X0 = x0*x0+b*x0+c;
arma::mat rs(n),x(n),xones(n),X(n),dx(n),out(n);
rs=pow((4*PI/3*n),-1./3.); //%# Added internal conversion to rs
x=sqrt(rs);
xones=arma::ones<arma::mat>(x.n_rows,x.n_cols);
X=x%x+b*x+xones*c;
dx=0.5/x; //%# Chain rule needs dx/drho!
out=dx%(2*X1/(rs%x)+A*(2./x-(2*x+b*xones)/X-4*b*xones/(Q*Q+(2.*x+b*xones)%(2.*x+b*xones))-(b*x0)/X0*(2/(x-x0*xones)-(2*x+b*xones)/X-4.*(2.*x0+b)/(Q*Q+(2*x+b*xones)%(2*x+b*xones)))));
out=(-rs/(3.*n))%out;
return out;

}

#endif // EXCVWN_H_INCLUDED
