#ifndef cIdag_H
#define cIdag_H

//! Hermitian transpose of backward fourier transformation operator
/*!
The hermitian transpose of fourier transformation from wavenumber space
to real space
 */

//#include "checkOperatorSize.h"
//#include "main.h"
#include <armadillo>
#include "latexComment.h"
#include "structs.h"

class paramStruct;


class cIdag
    {
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream

		//! The constructor
		/*!
		 \param Sinp a vector of dimensions of simulation domain
		 \param Numwavfunc ... number of wavefunctions to use
		 \param ltX ... pointer to latex class for documentation purposes */
        cIdag(const paramStruct Pa,arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX);//,checkOperatorSize<T> *chkPointer);
        ~cIdag();
		//! Laplace operator for input of values
		/*!
		 \param input as fourier coefficIdagents of wave functions expanded in plane waves
		 */
        arma::cx_mat operator*(arma::cx_mat input);
        arma::cx_mat operator*(arma::cx_mat* input);
    private:
        arma::Col<double> _S;
        arma::cx_mat _outputM;
        const paramStruct *_Pa;
        uvec _uv;

};


#endif // L_H
