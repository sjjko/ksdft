#include "sdOptimizer.h"

int sdOptimizer::setup(std::shared_ptr<arma::cx_mat> Wi)//,checkOperatorSize<T> *chkPointer)
{

    //! \brief orthogonalize the input wavefunction matrix
    //!
    //! \param Wi: input wavefunction matrix, complex, prod(S)*number_of_wavefunctions components
    //! computes $W \sqrt(W^{t}(O(W)))^{-1}$, where the square root acts component wise and $W^{t}(O(W))$
    //! is the integral value of each $\(\alpha,\beta\)$ component of the wavefunction
    std::string myFunctionString=" orthogonalize the wavefunction $W_{\\perp}=W \\sqrt(W^{t}(O(W)))^{-1}$ \\newline ";
    this->orthogonalizeWfunc(Wi); //orthogonalize W and return in W
    myLatexClass->setMyFunctionString(myFunctionString);

    return 1;
}

double sdOptimizer::optimize(std::shared_ptr<arma::cx_mat> Wi)
{
    verbosity(_Param,"optimize: we do an optimization for "+std::to_string(this->_sdNit)+" steps!",2,__FILE__,__LINE__);
    return this->optimize(this->_sdNit,Wi);

}

double sdOptimizer::optimize(const int Niterations,std::shared_ptr<arma::cx_mat> Wi)
{

    //! \brief steepest descent optimize function
    //!
    //! \param Niterations: number of iterations to take
    //! \param Wi: input wavefunction, complex, prod(S)*number_of_wavefunctions components

    stringstream myFunctionStringStream;
    myFunctionStringStream.precision(2);
    myFunctionStringStream<<" start of the steepest descent method \\newline " << endl;

    verbosity(_Param,"optimize Niterations: we do an optimization for "+std::to_string(Niterations)+" steps!",2,__FILE__,__LINE__);
    double E;
    //first we orthogonalize
    verbosity(_Param,"first set up the wavefunction",2,__FILE__,__LINE__);

    myFunctionStringStream<<" start of the steepest descent method \\newline " << endl;
    myFunctionStringStream<<" orthogonalize the wavefunction $W_{m}$ \\newline " << endl;
    myFunctionStringStream<<" $W=W \sqrt(W^{t}(O(W)))^{-1}$ \\newline " << endl;

    this->setup(Wi);
    verbosity(_Param,"calculate the initial electronic configuration energy of the system",2,__FILE__,__LINE__);

    myFunctionStringStream<<" compute the current energy of the system  \\newline " << endl;

    myFunctionStringStream<<getELatex(); //!< get the latex description of the getE energy computation function

    double E0=getE(this->_Op,this->_Param,*Wi.get(),this->_Vdual);
    int n=1;
    verbosity(_Param,"now do the sd optimization loop!",2,__FILE__,__LINE__);
    while(n<Niterations)
        {
            if(n==1) myFunctionStringStream<<getgradLatex(); //!< get the latex description of the getgrad energy gradient computation function
            arma::cx_mat g0=getgrad(this->_Op,this->_Param,*Wi.get(),this->_Vdual);
            if(n==1) myFunctionStringStream<<" now update the wavefunction by descending by alpha = " << this->_Param.alpha << " \\newline " << endl;
            if(n==1) myFunctionStringStream<<" $W=W-\\alpha \\nabla_{W} E$ \\newline " << endl;
            *Wi=*Wi-this->_alpha*g0;
            if(n==1) myFunctionStringStream<<" recalculate the energy \\newline " << endl;
            E=getE(this->_Op,this->_Param,*Wi,this->_Vdual);
            cout << "" << endl;
            cout << "====================================" << endl;
            cout << "sd: step " << n << " of " << Niterations << endl;
            cout << "sd: energy is:" << E << endl;
            cout << "sd: energy change is:" << (E-E0)/E0 << endl;
            cout << "====================================" << endl;

            if(n==1) myFunctionStringStream<<" \\subsubsection{  start the sd descent } " << endl;
            if(n==(Niterations-1))
            {
                myFunctionStringStream<<" after" << n << " steps " << endl;
                myFunctionStringStream<<" we have computed an energy of " << E << " hartree!" << endl;
            }

            E0=E;
            n++;


        }
    myLatexClass->setMyFunctionString(myFunctionStringStream.str());

    cout << "sd: final energy after " << Niterations << " iterations is:" << E << endl;
    return E;

}

sdOptimizer::sdOptimizer(const operatorStruct Operator,const paramStruct Parameter,const arma::mat Vdual,latexComment *ltX)
:optimizeBase(Operator,Parameter,Vdual)
{
    //! \brief initialize the steepest descent method class
    //!
    //! \param Operator: operators used, in one struct
    //! \param Parameter: parameters used, in one struct
    //! \param Vdual: the ion potential matrix in real space
    //! \param ltX: global latex class

//    this->_Nit=Parameter.Nit; //!< set number of iterations for sd
    this->_alpha=Parameter.alpha; //!< set alpha step parameter
    cout << "sd: intialize: Nit " << _sdNit << " alpha: " << _alpha << endl;

    myLatexClass=ltX;

    std::string myFunctionString=" \\textbf{\\large{}steepest descent solution method class} \\newline ";

    myLatexClass->setMyFunctionString(myFunctionString);

}

sdOptimizer::~sdOptimizer()
{
}




