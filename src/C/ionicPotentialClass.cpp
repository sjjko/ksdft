#include "ionicPotentialClass.h"

using namespace arma;
using namespace myFunctions;

ionicPotentialClass::ionicPotentialClass(const operatorStruct Op,paramStruct Pa,mat dr,cx_mat Sfi,arma::mat Xi,arma::mat G2i,string caseName)
{
    _caseName=caseName;
    _Op=Op;
    _Pa=Pa;
    _dr=dr;
    _Sf=Sfi;
    _X=Xi;
    _G2=G2i;
    _ListOfPotentials.push_back("harmonicOscillator");
    _ListOfPotentials.push_back("Hatom");
    _ListOfPotentials.push_back("Hmolecule");
    _ListOfPotentials.push_back("Germanium");
    verbosity(_Pa,"ionicPotentialClass::ionicPotentialClass finished initialization",2,__FILE__,__LINE__);

}

ionicPotentialClass::~ionicPotentialClass()
{
}

void ionicPotentialClass::getListOfPotentialPotentials()
{
    cout << "===================================================== " << endl;
    cout << "The list of available ionic potentials is as follows: " << endl;
    for (std::vector<string>::iterator it = _ListOfPotentials.begin() ; it != _ListOfPotentials.end(); ++it)
        {
        cout << *it << endl;
        }
    cout << "===================================================== " << endl;

}

int ionicPotentialClass::computePotential()
{

myFunctions::verbosity(_Pa,"ionicPotentialClass::getPotential: get the ionic potential for case "+_caseName,2,__FILE__,__LINE__);

if(_caseName=="harmonicOscillator")
    {
    double omega=2.;
    arma::mat V=0.5*pow(omega,2)*(_dr%_dr);
    this->_Vdual=real(*_Op.Jd*(*_Op.O*(*_Op.J*V)));
    }
else if(_caseName=="Hatom")
    {
    arma::mat Vps(_Pa.prodS,1,fill::zeros);
    Vps.col(0)=(-4.*arma::datum::pi*_Pa.Z/_G2);
    Vps(0)=0;
    this->_Vdual=real(*_Op.J*(Vps%_Sf));
    }
else if(_caseName=="Germanium")
    {
    verbosity(_Pa,"ionicPotentialClass::computePotential setup germanium ionic potential",2,__FILE__,__LINE__);
    double Z=_Pa.Z;
    double lambda=18.5;
    double rc=1.052;
    mat Gm=sqrt(_G2);
    checkArbitraryMatrix(Gm,"::computePotential Gm");

    mat Vps=-2*PI*exp(-PI*Gm/lambda)%cos(Gm*rc)%(Gm/lambda)/(1.-exp(-2.*PI*Gm/lambda));
    checkArbitraryMatrix(exp(-PI*Gm/lambda),"::computePotential exp(-PI*Gm/lambda)");
    checkArbitraryMatrix(exp(-PI*Gm/lambda)%cos(Gm*rc)%(Gm/lambda),"::computePotential exp(-PI*Gm/lambda)%cos(Gm*rc)%(Gm/lambda)");
    cout << find_nonfinite(Vps) << endl;
    cout << Gm(find_nonfinite(Vps)) << endl;
    //checkArbitraryMatrix(1./(1.-exp(-2.*Gm)),"::computePotential 1./(1.-exp(-2.*Gm))");
    Vps(0)=0;

    checkArbitraryMatrix(Vps,"::computePotential Vps1");

    verbosity(_Pa,"ionicPotentialClass::computePotential do the sum",2,__FILE__,__LINE__);
    for(int n=0;n<5;n++)
    {
        Vps=Vps+pow(-1.,n)*exp(-lambda*rc*n)/(1+square(n*lambda/Gm));
    }
    //checkArbitraryMatrix(Vps,"::computePotential Vps2");
    Vps(0)=0;

    myFunctions::cassert(((Vps.is_finite()) && (!Vps.has_nan())),ISCRITICAL,"ionicPotentialClass::getPotential: found nan or infinite value in Vps 2",__FILE__,__LINE__);

    Vps=Vps%(4.*PI*Z/square(Gm)*(1+exp(-lambda*rc)))-4.*PI*Z/square(Gm);
    //checkArbitraryMatrix(Vps,"::computePotential Vps3");
    Vps(0)=0;

    vec n={1,2,3,4};
    vec powM1n={-1.,1,-1,1};
    Vps(1)=4.*PI*Z*(1+exp(-lambda*rc))*(pow(rc,2)/2.+1./pow(lambda,2)*(PI*PI/6.+sum(powM1n%exp(-lambda*rc*n)/pow(n,2))));
    myFunctions::cassert(((Vps.is_finite()) && (!Vps.has_nan())),ISCRITICAL,"ionicPotentialClass::getPotential: found nan or infinite value in Vps 3",__FILE__,__LINE__);
    checkArbitraryMatrix(Vps,"::computePotential Vps4");

    verbosity(_Pa,"ionicPotentialClass::computePotential now transform Ge potential back to real space",2,__FILE__,__LINE__);
    //Vps(0)=0;
    this->_Vdual=real(*_Op.J*(Vps%_Sf));
    }
else
    {
    myFunctions::cassert(1==0,ISCRITICAL,"ionicPotentialClass::getPotential: potential for case "+_caseName + " not yet implemented!",__FILE__,__LINE__);
    }
    return 0;

    verbosity(_Pa,"ionicPotentialClass::computePotential Let's check the potential before we go on!",2,__FILE__,__LINE__);
    this->checkPotential();

    verbosity(_Pa,"ionicPotentialClass::computePotential finished",2,__FILE__,__LINE__);


}

void ionicPotentialClass::checkPotential()
{
    myFunctions::cassert(((this->_Vdual.is_finite()) && (!this->_Vdual.has_nan())),ISCRITICAL,"ionicPotentialClass::getPotential: found nan or infinite value in dual potential",__FILE__,__LINE__);
}

void ionicPotentialClass::checkArbitraryMatrix(mat M,string errorMessage)
{
    cout << " check matrix, got following indices for nan values: " << find_nonfinite(M) << endl;
    myFunctions::cassert(((M.is_finite()) && (!M.has_nan())),ISCRITICAL,"ionicPotentialClass: found nan or infinite value in "+errorMessage,__FILE__,__LINE__);
}

arma::mat ionicPotentialClass::getPotential()
{
return this->_Vdual;
}

