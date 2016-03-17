#ifndef ATOMICSYSTEM_H
#define ATOMICSYSTEM_H

#include "main.h"
#include "customAssert.h"
#include "pccgOptimizer.h"
#include "sdOptimizer.h"
#include "structs.h"
#include "inputParser.h"
#include <dlfcn.h>
#include "gnuPlotPlotting.h"
#include "ionicPotentialClass.h"

//! master class for setup and performance of dft computations

using namespace std;
using namespace arma;

class atomicSystem
{

//! \brief class which contains all ingredients to solve nonlinear energy functional from parser to postprocessing routines
//! is a wrapper for all the routines amassed so far

    public:

        atomicSystem(operatorStruct Op,paramStruct Pa,gnuPlotPlotting *gpLT,latexComment *latX);
        virtual ~atomicSystem();
        int getAtomicCoordinates();
        int setupGeometry();
        int setupPotential();
        int setupWavefunction(); //!<setup random wavefunction matrix as initial condition
        int postProcessing();
        int setupOptimizers();
        void plotX();
        string solveIt();
        int postVdual(const string texCaption);
        inline int orthogonalizeWavefunction()
            {
            //! \brief orthogonalize the wavefunction
            //! \param W - wavefunction as prodSxNo wavefunction complex matrix
            *_W=*_W*inv(sqrt((_W->t())*(*_Op.O*(*_W))));
            myFunctions::cassert(((this->_W->is_finite()) && (!this->_W->has_nan())),ISCRITICAL,"atomicSystem::orthogonalizeWavefunction: found nan or infinite value in wavefunction",__FILE__,__LINE__);
            return 0;
            }
        string postPsi();
        int postDensity(const string texCaption,const string fileNameEnding);
        inline string doPostProcessing(const string texCaption,const string fileNameEnding)
        {
            int returnInt=0;

            string returnString=this->postPsi();
            returnInt=this->postDensity(texCaption,fileNameEnding);

            return returnString;
        }


    private:

        string _caseName; //!< name of current atomic system
        paramStruct _Pa; //!< pointer to parameter class
        operatorStruct _Op; //!< pointer to operator class
        mat _mat_center_of_cell; //!< center of cell coordinates - convieniently (wastefully) stored in  prodSx3 matrix
        mat _r; //!< coordinate of computational cell: prodSx3 matrix
        mat _dr; //!< distance of computational cell from center: prodSx1 matrix
        mat _G; //!< wavenumber: prodSx3 matrix
        mat _G2; //!< squared wavenumber: prodSx1 matrix
        mat _Xt; //!< the atomic coordinate matrix
        latexComment *_ltX; //!< pointer to global latex documentation class
        mat _X; //!< atomic coordinates matrix transposed
        arma::cx_mat _Sf; //!< the structure function as exp(i G X), an prod(S)x1 matrix
        double _ionSelfEnergyUewald; //! self energy of ion cores

        mat _Vdual;

        std::shared_ptr<cx_mat> _W; //!< the wavefunction stored in complex prodSx3 matrix

        std::shared_ptr<pccgOptimizer> _pccg;
        std::shared_ptr<sdOptimizer> _sd;

        gnuPlotPlotting *_gpLT;

};

#endif // ATOMICSYSTEM_H
