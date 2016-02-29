#include "pccgOptimizer.h"

int pccgOptimizer::setup(std::shared_ptr<arma::cx_mat> Wi)//,checkOperatorSize<T> *chkPointer)
{
    return this->orthogonalizeWfunc(Wi); //orthogonalize W and return in W
}

//preconditioning
arma::cx_mat pccgOptimizer::K(const arma::cx_mat Wi)
{
    arma::cx_mat Wout(Wi);
    arma::mat ones(_G2);
    ones.fill(1.);
    for(int i=0;i<Wout.n_cols;i++)
    {
        Wout.col(i)=(ones/(ones+this->_G2)%(Wout.col(i)));
    }
    return Wout;
}



double pccgOptimizer::optimize(const int Niterations,std::shared_ptr<arma::cx_mat> Wi)
{
    double Energy_at_this_level;
    //first we orthogonalize
    verbosity(_Param,"first do the setup",2,__FILE__,__LINE__);
    this->setup(Wi);
    verbosity(this->_Param,"calculate the inital energy E0",2,__FILE__,__LINE__);
    double E0=getE(this->_Op,this->_Param,*Wi,this->_Vdual);
    verbosity(this->_Param,"compute the initial energy gradient",2,__FILE__,__LINE__);
    arma::cx_mat g_nm1=getgrad(_Op,_Param,*Wi,_Vdual);
    verbosity(this->_Param,"do the preconditioning",2,__FILE__,__LINE__);
    arma::cx_mat d_nm1=-this->K(g_nm1);
    verbosity(this->_Param,"now start looping",2,__FILE__,__LINE__);

    int n=0;
    while(n<Niterations)
        {
            arma::cx_mat g_n=getgrad(this->_Op,this->_Param,*Wi,this->_Vdual);
            if(n>1)
            {
                double acosine=Prod(g_n,d_nm1)/sqrt((Prod(g_n,g_n)*Prod(d_nm1,d_nm1)));
                double cgtest=Prod(g_n,this->K(g_nm1))/sqrt((Prod(g_n,this->K(g_n))*Prod(g_nm1,this->K(g_nm1))));
            }
            double beta=Prod(g_n,this->K(g_n))/Prod(g_nm1,this->K(g_nm1));
            arma::cx_mat d_n=-K(g_n)+beta*d_nm1;
            arma::cx_mat Wtrial=*Wi+this->_alpha*d_n;
            arma::cx_mat gtrial=getgrad(this->_Op,this->_Param,Wtrial,this->_Vdual);
            double alpha=this->_alpha*Prod(g_n,d_n)/Prod(g_n-gtrial,d_n);
            *Wi=*Wi+alpha*d_n;
            Energy_at_this_level=getE(this->_Op,this->_Param,*Wi,this->_Vdual);
            g_nm1=g_n;
            d_nm1=d_n;
            cout << "====================================" << endl;
            cout << "pccg: step " << n << " of " << Niterations << endl;
            cout << "pccg: energy is:" << Energy_at_this_level << endl;
            cout << "pccg: energy change is:" << (Energy_at_this_level-E0)/E0 << endl;
            cout << "====================================" << endl;
            E0=Energy_at_this_level;
            n++;
        }

    cout << "pccg: final energy after " << Niterations << " iterations is:" << Energy_at_this_level << endl;
    return Energy_at_this_level;

}

double pccgOptimizer::optimize(std::shared_ptr<arma::cx_mat> Wi)
{

double E=this->optimize(this->_pccgNit,Wi);
return E;

}

pccgOptimizer::pccgOptimizer(const operatorStruct Operator,const paramStruct Parameter,const arma::mat Vdual,const arma::mat G2Input,latexComment *ltX)
:optimizeBase(Operator,Parameter,Vdual)
{
    //this->_Nit=Parameter._PccgNit;
    this->_alpha=Parameter.alpha;
    if(_alpha==0){_alpha=0.5e-3;}
    _G2=G2Input;

    myLatexClass=ltX;

    std::string myFunctionString=" \\textbf{\\large{}conjugate optimizer class} \\newline ";

    myLatexClass->setMyFunctionString(myFunctionString);
}

/*pccgOptimizer::~pccgOptimizer()
{
}*/




