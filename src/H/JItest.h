#ifndef JITEST_H_INCLUDED
#define JITEST_H_INCLUDED


#include "main.h"

bool IJtest(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid);
bool IJtestMultiColumn(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid);

bool IJtest(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid)
{
    //! \brief do a J I inversion transformation test
    //!
    //! \param O an operator struct
    //! \param P a parameter struct
    //! \return S, the S column vector containing the grid dimension

bool returnval=false;
arma::cx_mat dW=arma::randu<arma::cx_mat>(arma::prod(Sgrid),1);
//cout << dW << endl;
arma::cx_mat result(*O.J*(*O.I*dW)/dW);
//cout << real(sum(*O.J*(*O.I*dW)-dW));
double realdifference=as_scalar(real(sum(*O.J*(*O.I*dW)-dW)));
double imagdifference=as_scalar(imag(sum(*O.J*(*O.I*dW)-dW)));

if((realdifference+imagdifference)<0.01){returnval=true;}
return returnval;
}

bool IJtestMultiColumn(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid)
{
    //! \brief do a J I inversion transformation test with a multicolumn matrix input
    //!
    //! \param O an operator struct
    //! \param P a parameter struct
    //! \return S, the S column vector containing the grid dimension

bool returnval=false;
arma::cx_mat dW=arma::randu<arma::cx_mat>(arma::prod(Sgrid),5);
//cout << dW << endl;
arma::cx_mat result(*O.J*(*O.I*dW)/dW);
//cout << real(sum(*O.J*(*O.I*dW)-dW));
arma::vec rowSums=vectorise(real(sum(*O.J*(*O.I*dW)-dW)));
double realdifference=as_scalar(sum(rowSums,0));
rowSums=vectorise(imag(sum(*O.J*(*O.I*dW)-dW)));
double imagdifference=as_scalar(sum(rowSums,0));

if((realdifference+imagdifference)<0.01){returnval=true;}
return returnval;
}


#endif
