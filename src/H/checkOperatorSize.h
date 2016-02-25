#ifndef CHECKOPERATORSIZE_H
#define CHECKOPERATORSIZE_H

#include "../../main.H"

// checkOp class
/** check Operator class has a function which checks the size of two matrices operated on*/

template<class T> class checkOperatorSize
    {
    public:
        checkOperatorSize(){};
        ~checkOperatorSize();
        inline int checkOps(T op1,T op2)
        {
            try{
                std::cout << arma::size(op1) << " " << arma::size(op2) << std::endl;
                return 1;}
            catch(std::exception &e)
            {return -1;}

        };

};

template class checkOperatorSize<arma::vec>;


#endif // L_H
