#ifndef OPTIMIZERBASE_H
#define OPTIMIZERBASE_H

//! The base class for various optimizers implemented
/*!
 This class provides the basis for all the optimization algorithms implemented
 */

//#include "checkOperatorSize.h"
#include "main.h"
#include "structs.h"
#include "dft_functions.h"
#include "customAssert.h"
#include "structs.h"

class optimizeBase
    {
    public:
		//! The constructor
		/*!
		 \param G2inp a vector of wave vector squared
		 \param Rinp a 3x3 matrix of system size
		 \param chkPointer a pointer to a class checking the validity of operations */
        inline optimizeBase(const operatorStruct Operator,const paramStruct Parameter,const arma::mat Vdual)
                {
                this->_Op=Operator;
                this->_Param=Parameter;
                this->_Vdual=Vdual;
                _sdNit=Parameter.sdNit;
                _pccgNit=Parameter.pccgNit;
                };
        //~optimizeBase();
		//! Laplace operator for input of values
		/*!
		 \param input as fourier coefficients of wave functions expanded in plane waves
		 */
        inline virtual int setup() {return 1;}
        inline int orthogonalizeWfunc(std::shared_ptr<arma::cx_mat> Winp)
        {
            try
            {
                arma::cx_mat tempM=sqrt(arma::inv(Winp->t()*(*_Op.O*(Winp.get()))));
                arma::cx_mat inpM=*Winp;
                verbosity(this->_Param,"W+(O(W)) has nrows:"+std::to_string(tempM.n_rows)+" and ncols: "+std::to_string(tempM.n_cols),2,__FILE__,__LINE__);
                *Winp=inpM*tempM; //sqrt(arma::inv(Winp->t()*(*_Op.O*(*Winp))));
                }
                catch(exception& e)
                {
                cout << "error in optimizeBase in orthogonalizeWfunc!" << endl;
                cout << "Standard exception " << e.what() << "raised!" << endl;
                throw e;
                  }
                return 1;
        }
        /*inline optimization of given wavefunction (overwritten)*/
        inline virtual double optimize(std::shared_ptr<arma::cx_mat> W){return 1.0;}
        inline virtual double optimize(const int Niterations,std::shared_ptr<arma::cx_mat> Wi){return 1.0;};

        operatorStruct _Op;
        paramStruct _Param;
        arma::mat _Vdual;

    protected:
        double _alpha;
        int _sdNit;
        int _pccgNit;


};

#endif // L_H
