#ifndef TESTINGCLASS_H
#define TESTINGCLASS_H

#include "customAssert.h"
//#include "slice.h"
//#include "main.h"
#include "structs.h"
#include "latexComment.h"
#include "plotting.h"
#include "openCVPlotting.h"
#include "ppm.h"
#include "gnuPlotPlotting.h"
#include "dft_functions.h"
#include "sdOptimizer.h"
#include "pccgOptimizer.h"

class testingClass
{

    private:
        latexComment *_myLatexPtr;
        paramStruct _myP;

    public:
        testingClass(latexComment *latX,paramStruct P);
        virtual ~testingClass();
        bool IJtest(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid);
        bool IJtestMultiColumn(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid);
        bool poissonEquationTest(const operatorStruct O,const paramStruct P,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell);
        bool ewaldTest(const operatorStruct O,const paramStruct P,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell,const arma::cx_mat Sf,arma::mat X);
        bool checkHermitianConsistency(const operatorStruct O,const paramStruct P,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell,const arma::cx_mat Sf,arma::mat X);
        bool multicolumnTest(const operatorStruct O,const paramStruct P,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell,const arma::cx_mat Sf,arma::mat X);
        bool schroedingerTest(const operatorStruct Op,const paramStruct Pa,arma::Col<double> S,Mat<double> R,
        mat r,mat mat_center_of_cell,const arma::cx_mat Sf,arma::mat X,arma::Mat<double> G2);

        inline void testForIJ(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid)
        {if(this->IJtest(O,P,Sgrid)){cout << "J I inversion test passed" << endl;}else{cout << "J I inversion test failed" << endl;}};
        inline void testForIJMultiColumn(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid)
        {if(this->IJtestMultiColumn(O,P,Sgrid)){cout << "J I multicolumn inversion test passed" << endl;}else{cout << "J I multicolumn inversion test failed" << endl;}};
        inline void testPoisson(const operatorStruct O,const paramStruct P,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell)
        {
        if(this->poissonEquationTest(O,P,S,R,r,mat_center_of_cell))
            {
                cout << "Poisson test passed" << endl;
            }
            else
            {
                cout << "Poisson test failed" << endl;
            }
        };
        inline void testEwald(const operatorStruct O,const paramStruct P,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell,arma::cx_mat Sf,arma::mat X)
        {
        if(this->ewaldTest(O,P,S,R,r,mat_center_of_cell,Sf,X))
            {
                cout << "Ewald test passed" << endl;
            }
            else
            {
                cout << "Ewald test failed" << endl;
            }
        };
        inline void testSchroedinger(const operatorStruct Op,const paramStruct Pa,arma::Col<double> S,Mat<double> R,
        mat r,mat mat_center_of_cell,const arma::cx_mat Sf,arma::mat X,arma::Mat<double> G2)
        {
        if(this->schroedingerTest(Op,Pa,S,R,r,mat_center_of_cell,Sf,X,G2))
            {
                cout << "Schroedinger equation test passed" << endl;
            }
            else
            {
                cout << "Schroedinger equation test failed" << endl;
            }
        };
        inline void testHermitian(const operatorStruct O,const paramStruct P,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell,arma::cx_mat Sf,arma::mat X)
        {
        if(this->checkHermitianConsistency(O,P,S,R,r,mat_center_of_cell,Sf,X))
            {
                cout << "Hermitian test passed" << endl;
            }
            else
            {
                cout << "Hermitian test failed" << endl;
            }
        };
        inline void testMultiColumn(const operatorStruct O,const paramStruct P,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell,arma::cx_mat Sf,arma::mat X)
        {
        if(this->multicolumnTest(O,P,S,R,r,mat_center_of_cell,Sf,X))
            {
                cout << "multicolumn test passed" << endl;
            }
            else
            {
                cout << "multicolumn test failed" << endl;
            }
        };
        inline void testDiagOuter()
        {
        bool testPassed=false;
        cx_mat a(15,3);
        cx_mat b(15,3);
        a.set_real(randn<mat>(15,3));
        b.set_real(randn<mat>(15,3));
        a.set_imag(randn<mat>(15,3));
        b.set_imag(randn<mat>(15,3));

        //! check for maximal element in the two approaches
        std::complex<double> maxV=as_scalar(max(diagvec(a*b.t())-diagOuter(_myP,a,b)));
        const std::complex<double> bound(1.0e-5,1.0e-05);
        if(std::real(maxV*std::conj(maxV))<std::real(bound*std::conj(bound))){testPassed=true;}

        if(testPassed)
            {
                cout << "diagOuter test passed" << endl;
            }
            else
            {
                cout << "diagOuter test failed" << endl;
            }

        };
        inline void testOrthogonalizationOfWF(const operatorStruct Op,const paramStruct Pa,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell,arma::cx_mat Sf,arma::mat X)
        {
        //! \brief test the orthogonalization procedure for Wavefunctions


            arma_rng::set_seed(1.0); //!< not so random as it seems
            arma::cx_mat Wtrial=arma::randu<arma::cx_mat>(prod(S),Pa.number_of_wavefunctions);

            //Wtrial.set_real(arma::randu<arma::mat>(prodS,Pa.number_of_wavefunctions));
            //Wtrial.set_imag(arma::randu<arma::mat>(prodS,Pa.number_of_wavefunctions));
            //cout << Wtrial << endl;
            //cout << size(Wtrial) << endl;
            arma::cx_mat squareM=Wtrial.t()*(*Op.O*Wtrial);
            arma::cx_mat Y=Wtrial*arma::sqrtmat(arma::inv(squareM));
            //cout << Wtrial.row(0) << endl;
            arma::cx_mat OrthoMatrix=trans(Y)*(*Op.O*Y);

            arma::cx_mat Ui=Uinvers(Op,Pa,Wtrial);
            Y=Wtrial*sqrtmat(Ui);
            arma::cx_mat OrthoMatrix2=trans(Y)*(*Op.O*Y);

            cout << OrthoMatrix2 - OrthoMatrix << endl;
            cout << OrthoMatrix2 << endl;


        }


};

#endif // TESTINGCLASS_H
