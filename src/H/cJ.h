#ifndef CJ_H
#define CJ_H

//! forward fourier transformation operator class
/*!
 This class provides operators for computing the transformation
 of quantities from real to fourier space
 */

#include <armadillo>
#include "latexComment.h"

class cJ
    {
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream

        cJ(arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX);
		//! The constructor
        /*!
         \param Sinp a vector of dimensions of simulation domain
		 \param Numwavfunc ... number of wavefunctions to use
		 \param ltX ... pointer to latex class for documentation purposes !*/
        ~cJ();
        arma::cx_mat operator*(arma::mat input); //! fourier transformation of real matrix
        arma::cx_mat operator*(arma::cx_mat input); //! fourier transformation class of complex matrix
        arma::cx_mat operator*(arma::cx_mat* input); //! fourier transformation class of complex matrix represented by pointer
    private:
        arma::Col<double> _S;
        arma::cx_mat _outputM;
};

#endif // L_H
