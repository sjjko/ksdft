#ifndef LINV_H
#define LINV_H

#include <armadillo>
#include "latexComment.h"

///inverse Laplacian
/** Linv: compute the Laplacian of a number N of wavefunction coefficients */

class Linv
//: public checkOperatorSize<T>
{
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream
        Linv(const arma::mat G2inp,const arma::mat Rinp,const int Numwavfunc,latexComment *ltX); //constructor: input G2 squared wavevector array, R matrix contains size of computing box
        ~Linv();

//        arma::cx_vec  operator*(arma::cx_vec input);
		arma::cx_mat  operator*(arma::cx_mat input); //overloaded op. for copied input
        arma::cx_mat  operator*(arma::cx_mat* input); // ovrl. op. for refer. input
    private:
        arma::mat _G2i;
        arma::mat _Faktor;
        arma::mat _Ri;
        arma::cx_mat  _outputM;
		//checkOperatorSize<T> *_chk;
		//        std::auto_ptr<arma::mat> Rint;
};



#endif // L_H
