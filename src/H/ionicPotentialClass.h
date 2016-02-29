#ifndef IONICPOTENTIALCLASS_H
#define IONICPOTENTIALCLASS_H

#include "main.h"
#include "structs.h"
#include "latexComment.h"
#include "customAssert.h"
//#include <armadillo>

using namespace arma;

class ionicPotentialClass
{
    public:
        ionicPotentialClass(operatorStruct Op,paramStruct Pa,mat dr,cx_mat Sfi,arma::mat Xi,arma::mat G2i,string caseName);
        virtual ~ionicPotentialClass();
        arma::mat getPotential(); //! < retrieve dual potential for case specified in string
        void checkPotential();
        int computePotential();
        void getListOfPotentialPotentials();
        void checkArbitraryMatrix(mat M,string errorMessage);

    private:
        arma::mat _Vdual; //! the potential to return for further computations

        string _caseName;
        operatorStruct _Op;
        paramStruct _Pa;
        mat _dr;
        cx_mat _Sf;
        arma::mat _X;
        arma::mat _G2;
        std::vector<std::string> _ListOfPotentials;

};

#endif // IONICPOTENTIALCLASS_H
