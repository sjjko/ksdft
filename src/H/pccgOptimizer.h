#ifndef PCCGOPTIMIZER_H
#define PCCGOPTIMIZER_H

//! pccg optimizer implementation
/*!
 conjugate gradient optimization algorithm
 */

//#include "checkOperatorSize.h"
#include "main.h"
#include "optimizerBase.h"
#include "dft_functions.h"
#include "latexComment.h"

class pccgOptimizer:public optimizeBase
    {
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream
        pccgOptimizer(const operatorStruct Operator,const paramStruct Parameter,const arma::vec Vreal,const arma::mat G2Input,latexComment *ltX);
//        ~pccgOptimizer();
		//! Laplace operator for input of values
		/*!
		 \param pointer to wavefunction coefficient matrix (is altered by setup)
		 */
        int setup(std::shared_ptr<arma::cx_mat> W);
        arma::cx_mat K(const arma::cx_mat W);
        double optimize(const int  Niterations,std::shared_ptr<arma::cx_mat> W);
        double optimize(std::shared_ptr<arma::cx_mat> W);

    private:
        double _alpha;
        arma::mat _G2;
};

#endif // L_H
