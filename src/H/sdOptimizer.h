#ifndef SDOPTIMIZER_H
#define SDOPTIMIZER_H

//! The base class for various optimizers implemented
/*!
 This class provides the basis for all the optimization algorithms implemented
 */

//#include "checkOperatorSize.h"
#include "main.h"
#include "optimizerBase.h"
#include "dft_functions.h"
#include "latexComment.h"

class sdOptimizer:public optimizeBase
    {
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream
        sdOptimizer(const operatorStruct Operator,const paramStruct Parameter,const arma::mat Vreal,latexComment *ltX);
        ~sdOptimizer();
		//! Laplace operator for input of values
		/*!
		 \param pointer to wavefunction coefficient matrix (is altered by setup)
		 */
        int setup(std::shared_ptr<arma::cx_mat> W);
        /* optimization of given wavefunction (overwritten)*/
        double optimize(std::shared_ptr<arma::cx_mat> W);
        double optimize(const int Niterations,std::shared_ptr<arma::cx_mat> Wi);
        //virtual int orthogonalizeWfunc(std::shared_ptr<arma::cx_mat> Winp);
        //inline int orthogonalizeWfunc(std::shared_ptr<arma::cx_mat> Winp){};

//    private:
//        const double _alpha=0.5e-3;

};

#endif // L_H
