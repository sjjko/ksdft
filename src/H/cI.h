#ifndef CI_H
#define CI_H

//! backward fourier transformation operator class
/*!
 This class provides operators for computing the transformation
 of quantities from fourier to real space
 */

 #include "main.h"
#include <armadillo>
#include "latexComment.h"
#include "structs.h"
#include "customAssert.h"

class paramStruct;

class cI
    {
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream
        paramStruct *_Pa;

        cI(paramStruct Pa, arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX);//,checkOperatorSize<T> *chkPointer);
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
        uvec _uv;
        cx_mat _full;
};


#endif
