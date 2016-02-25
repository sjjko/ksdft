#ifndef CJ_H
#define CJ_H

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


class cJ
    {
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream
		//! The constructor
		/*!
		 \param G2inp a vector of wave vector squared
		 \param Rinp a 3x3 matrix of system size
		 \param chkPointer a pointer to a class checking the validity of operations */
        cJ(arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX);//,checkOperatorSize<T> *chkPointer);
        ~cJ();
		//! Laplace operator for input of values
		/*!
		 \param input as fourier coefficients of wave functions expanded in plane waves
		 */
        arma::cx_mat operator*(arma::mat input);
        arma::cx_mat operator*(arma::cx_mat input);
        arma::cx_mat operator*(arma::cx_mat* input);
        //arma::cx_vec operator*(arma::cx_vec input);
        //arma::cx_vec operator*(arma::vec input);
    private:
        arma::Col<double> _S;
        //checkOperatorSize<T> *_checkPointer;
        arma::cx_mat _outputM;

//        std::auto_ptr<arma::mat> Rint;

};


#endif // L_H
