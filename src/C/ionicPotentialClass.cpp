#include "ionicPotentialClass.h"

using namespace arma;
using namespace myFunctions;


double ionicPotentialClass::computeEwaldIonSelfEnergy()
{
//! \brief this function computes the self energy of the ions in their electrostatic field
double selfEnergy;


const double sigma1=0.25;

_myLatexPtr->newLine(" charge distribution 1 as: $g_{1}=Z{\\exp(-{dr^{2} \\over 2 \sigma_{1})^{2}}) \\over 2 \\pi \\sigma_{1}^{3}}$");

verbosity(_Pa,"ionicPotentialClass::computeEwaldIonSelfEnergy compute g1",2,__FILE__,__LINE__);

arma::cx_mat g1(_dr.n_rows,1,arma::fill::zeros);
g1.set_real(_Pa.Z*arma::exp(-(_dr%_dr)/(2.*pow(sigma1,2)))/sqrt(pow(2.*arma::datum::pi*sigma1*sigma1,3)));

verbosity(_Pa,"sum dr = "+std::to_string(arma::as_scalar(arma::sum(arma::sum(_dr)))),2,__FILE__,__LINE__);
verbosity(_Pa,"sum g1 = "+std::to_string(arma::as_scalar(real(arma::sum(arma::sum(g1))))),2,__FILE__,__LINE__);

cassert(!_dr.has_inf(),ISCRITICAL,"dr has inf",__FILE__,__LINE__);
cassert(!_dr.has_nan(),ISCRITICAL,"dr has nan",__FILE__,__LINE__);

cassert(!g1.has_inf(),ISCRITICAL,"g1 has inf",__FILE__,__LINE__);
cassert(!g1.has_nan(),ISCRITICAL,"g1 has nan",__FILE__,__LINE__);


//arma::cx_mat n(g1.n_elem,1);
//n.fill(0);
verbosity(_Pa,"ionicPotentialClass::computeEwaldIonSelfEnergy: now extract density by applying the charge distribution to all the fields using the structure function Sf",2,__FILE__,__LINE__);

_myLatexPtr->newLine(" charge density n is $n=g_{2}-g_{1}$");

arma::cx_mat nc(arma::size(_dr),fill::zeros);

cassert(!_Sf.has_inf(),ISCRITICAL,"Sf has inf",__FILE__,__LINE__);
cassert(!_Sf.has_nan(),ISCRITICAL,"Sf has nan",__FILE__,__LINE__);
cassert(_Sf.is_finite(),ISCRITICAL,"Sf has nan",__FILE__,__LINE__);

nc.set_real(real(*_Op.I*((*_Op.J*g1)%_Sf)));

cassert(!nc.has_inf(),ISCRITICAL,"nc has inf",__FILE__,__LINE__);
cassert(!nc.has_nan(),ISCRITICAL,"nc has nan",__FILE__,__LINE__);
cassert(nc.is_finite(),ISCRITICAL,"nc has nan",__FILE__,__LINE__);

//gpL.plotMatrix3Slices("ewaldrealSf",real(SfM),S);
//gpL.plotMatrix3Slices("ewaldimagSf",imag(SfM),S);


arma::mat nreal=real(nc);

//gpL.plotMatrix3Slices("nreal",nreal,S);

verbosity(_Pa,"ewald: compute density distribution using structure function Sf",2,__FILE__,__LINE__);

verbosity(_Pa,"ionicPotentialClass::computeEwaldIonSelfEnergy: max nc = "+std::to_string(arma::as_scalar(real(arma::max(arma::max(nc))))),2,__FILE__,__LINE__);
verbosity(_Pa,"ionicPotentialClass::computeEwaldIonSelfEnergy: min nc = "+std::to_string(arma::as_scalar(real(arma::min(arma::min(nc))))),2,__FILE__,__LINE__);

verbosity(_Pa,"ionicPotentialClass::computeEwaldIonSelfEnergy: max nreal = "+std::to_string(arma::as_scalar(arma::max(arma::max(nreal)))),2,__FILE__,__LINE__);
verbosity(_Pa,"ionicPotentialClass::computeEwaldIonSelfEnergy: min nreal = "+std::to_string(arma::as_scalar(arma::min(arma::min(nreal)))),2,__FILE__,__LINE__);

verbosity(_Pa,"ionicPotentialClass::computeEwaldIonSelfEnergy:  nc is a real stored in a complex matrix",2,__FILE__,__LINE__);
verbosity(_Pa,"ionicPotentialClass::computeEwaldIonSelfEnergy:  now compute potential phi",2,__FILE__,__LINE__);

arma::cx_mat phik = solvePoisson(_Op,_Pa,nreal);
arma::cx_mat phi = *_Op.I*phik;

double Unum=0.5*real(arma::accu(((*_Op.J*phi).t()*(*_Op.O*(*_Op.J*nc)))));

arma::mat Xtranspose=_X.t();
verbosity(_Pa,"ionicPotentialClass::computeEwaldIonSelfEnergy:  we have "+std::to_string(Xtranspose.n_rows)+" atoms in our system!",2,__FILE__,__LINE__);

double Uself=_Pa.Z*_Pa.Z/(2.*sqrt(arma::datum::pi))*(1./sigma1)*Xtranspose.n_rows;


selfEnergy=Unum-Uself;

return selfEnergy;
}


ionicPotentialClass::ionicPotentialClass(const operatorStruct Op,paramStruct Pa,mat dr,mat mat_center_of_cell,cx_mat Sfi,arma::mat Xi,arma::mat G2i,string caseName,latexComment *latX)
{
    //! \brief setup of geometric and physical quantities for potential definition

    //! \param O ... an operator struct
    //! \param Pa ... a parameter struct
    //! \param dr ... NoGridCellx1 matrix holding distances to cell center
    //! \param Sfi ... complex NoGridCellx1 matrix holding structure factor Sfi(k)=sum_atom exp(i G X_atom)
    //! \param Xi ... real Natom matrix holding positions of the atomic cores
    //! \param G2i ... real NoGridCellx1 squared wavevector matrix
    //! \param caseName ... string describing name of the case at hand

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
    this->_myLatexPtr=latX;
    _mat_center_of_cell=mat_center_of_cell;

}

ionicPotentialClass::~ionicPotentialClass()
{
}

void ionicPotentialClass::getListOfPotentialPotentials()
{
    //! \brief output all potentials defined on the screen


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

    //! \brief computational module - here all the potentials in use are defined


myFunctions::verbosity(_Pa,"ionicPotentialClass::getPotential: get the ionic potential for case "+_caseName,2,__FILE__,__LINE__);

if(_caseName=="harmonicOscillator" || _caseName=="quantumDot")
    {
    double omega=2.;
    arma::mat V=0.5*pow(omega,2)*(_dr%_dr);
    this->_Vdual=real(*_Op.Jd*(*_Op.O*(*_Op.J*V)));
    }
else if(_caseName=="Hatom" || _caseName=="Hmolecule")
    {
    arma::mat Vps(_G2.n_rows,1,fill::zeros);
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

    //Vps=Vps.*4*pi*Z./Gm.^2*(1+exp(-lambda*rc))-4*pi*Z./Gm.^2;

    Vps=Vps*4.*PI*Z/square(Gm)*(1+exp(-lambda*rc))-4.*PI*Z/square(Gm);
    //checkArbitraryMatrix(Vps,"::computePotential Vps3");
    Vps(0)=0;

    vec n={1,2,3,4};
    vec powM1n={-1.,1,-1,1};
    Vps(0)=4.*PI*Z*(1+exp(-lambda*rc))*(pow(rc,2)/2.+1./pow(lambda,2)*(PI*PI/6.+sum(powM1n%exp(-lambda*rc*n)/pow(n,2))));
    myFunctions::cassert(((Vps.is_finite()) && (!Vps.has_nan())),ISCRITICAL,"ionicPotentialClass::getPotential: found nan or infinite value in Vps 3",__FILE__,__LINE__);
    checkArbitraryMatrix(Vps,"::computePotential Vps4");

    verbosity(_Pa,"ionicPotentialClass::computePotential now transform Ge potential back to real space",2,__FILE__,__LINE__);
    //Vps(0)=0;
    this->_Vdual=real(*_Op.J*(Vps%_Sf));

    gnuPlotPlotting gplt;
    string name=_Pa.caseName+"_Vdualdirekt";
    string imgName=gplt.plotAMatrixSlice(_Pa,name,this->_Vdual,_Pa.S,0);
    imgName=gplt.plotAMatrixSlice(_Pa,name,this->_Vdual,_Pa.S,1);
    imgName=gplt.plotAMatrixSlice(_Pa,name,this->_Vdual,_Pa.S,2);

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
    //! \brief check computed potential for finiteness and naNs


    myFunctions::cassert(((this->_Vdual.is_finite()) && (!this->_Vdual.has_nan())),ISCRITICAL,"ionicPotentialClass::getPotential: found nan or infinite value in dual potential",__FILE__,__LINE__);
}

void ionicPotentialClass::checkArbitraryMatrix(mat M,string errorMessage)
{
    //! \brief check a matrix for finiteness and naNs


    cout << " check matrix, got following indices for nan values: " << find_nonfinite(M) << endl;
    myFunctions::cassert(((M.is_finite()) && (!M.has_nan())),ISCRITICAL,"ionicPotentialClass: found nan or infinite value in "+errorMessage,__FILE__,__LINE__);
}

arma::mat ionicPotentialClass::getPotential()
{
    //! \brief return computed dual potential as real matrix

return this->_Vdual;
}

