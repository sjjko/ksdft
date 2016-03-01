#ifndef cJdag_H
#define cJdag_H

//! Hermitian transpose of forward fourier transformation operator
/*!
The hermitian transpose of fourier transformation from real to wavenumber space
 */

#include <armadillo>
#include "latexComment.h"

class cJdag
    {
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream

         /*!
         \param Sinp a vector of dimensions of simulation domain
		 \param Numwavfunc ... number of wavefunctions to use
		 \param ltX ... pointer to latex class for documentation purposes !*/

        cJdag(arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX); //! The constructor
        ~cJdag();

        arma::cx_mat operator*(arma::cx_mat input); //! hermitian fourier transformation operator matrix input
        arma::cx_mat operator*(arma::cx_mat* input); //! hermitian fourier transformation operator matrix pointer input
    private:
        arma::Col<double> _S; //! matrix carrying the dimension of the simulation domain
        arma::cx_mat _outputM;
};

#endif // L_H
