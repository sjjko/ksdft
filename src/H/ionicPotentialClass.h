#ifndef IONICPOTENTIALCLASS_H
#define IONICPOTENTIALCLASS_H

//! setup of ionic potentials - output is Vdual

#include "main.h"
#include "structs.h"
#include "latexComment.h"
#include "customAssert.h"
#include "gnuPlotPlotting.h"
#include "dft_functions.h"
//#include <armadillo>


using namespace arma;

class ionicPotentialClass
{
    public:
        ionicPotentialClass(operatorStruct Op,paramStruct Pa,mat dr,mat mat_center_of_cell,cx_mat Sfi,arma::mat Xi,arma::mat G2i,string caseName,latexComment *latX);
        //! \param Op .. operator structure holds
        virtual ~ionicPotentialClass();
        arma::mat getPotential(); //! < retrieve dual potential for case specified in string
        void checkPotential();
        int computePotential();
        void getListOfPotentialPotentials();
        void checkArbitraryMatrix(mat M,string errorMessage);
        double computeEwaldIonSelfEnergy();

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
        latexComment *_myLatexPtr; //!< pointer to global latex documentation class
        mat _mat_center_of_cell; //!< center of cell coordinates - convieniently (wastefully) stored in  prodSx3 matrix

};

#endif // IONICPOTENTIALCLASS_H
