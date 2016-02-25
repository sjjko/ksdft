#ifndef CI_H
#define CI_H

//! The Laplace operator class
/*!
 This class provides operators for computing the Laplacian
 It can be applied on fourier transformed sets of wave functions
 means on sets of fourier coefficients of plane wave functions
 */

//#include "checkOperatorSize.h"
//#include "main.h"
#include <armadillo>
#include "latexComment.h"

class cI
    {
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream

        cI(arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX);//,checkOperatorSize<T> *chkPointer);
		//! The constructor
		/*!
		 \param Sinp vector containing the number of latticepoints along each dimension
		 \param Numwavfunc number of wavefunctions
		 \param ltX pointer to latex class instance - for documentation */
        ~cI();
        arma::cx_mat operator*(arma::cx_mat input);
        arma::cx_mat operator*(arma::cx_mat* input);
    private:
        arma::Col<double> _S;
        //checkOperatorSize<T> *_checkPointer;
        arma::cx_mat _outputM;
};


#endif
