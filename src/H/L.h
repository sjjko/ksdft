#ifndef L_H
#define L_H

//! The Laplace operator class
/*!
 This class provides operators for computing the Laplacian
 It can be applied on fourier transformed sets of wave functions
 means on sets of fourier coefficients of plane wave functions
 */

//#include "checkOperatorSize.h"
#include "main.h"
#include "latexComment.h"

class Lop //: public latexComment
    {

    public:

        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream
        Lop(const arma::mat G2inp,const arma::mat Rinp,const int Numwavfunc,latexComment *ltX);//,checkOperatorSize<T> *chkPointer);
		//! The constructor
		/*!
		 \param G2inp a vector of wave vector squared
		 \param Rinp a 3x3 matrix of system size
		 \param Numwavfunc number of wavefunctions to compute
		 \param ltX - pointer to a latex comment class instance */

        ~Lop();
		//! Laplace operator for input of values
		/*!
		 \param input as fourier coefficients of wave functions expanded in plane waves
		 */
        arma::cx_mat operator*(arma::cx_mat input);
        arma::cx_mat operator*(arma::cx_mat* input);
    private:
        arma::mat _G2i;
        arma::mat _Ri;
        arma::mat _Faktor;
        arma::cx_mat _outputM;
        //checkOperatorSize<T> *_checkPointer;

//        std::auto_ptr<arma::mat> Rint;

};


#endif // L_H
