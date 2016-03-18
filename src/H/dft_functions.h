#ifndef DFTFUNCTIONS_H_INCLUDED
#define DFTFUNCTIONS_H_INCLUDED

#include "main.h"
#include "excVWN.h"
#include "excpVWN.h"
#include "customAssert.h"
#include "structs.h"

using namespace arma;


double getE(const operatorStruct O,const paramStruct P,arma::cx_mat Wi,arma::mat VdualIon);
arma::cx_mat getgrad(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi,const arma::mat VdualIon);
arma::cx_mat Uinvers(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi);
arma::mat computeDensityFromWavefuncs(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi,const arma::cx_mat Ui);
arma::cx_mat HW(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi,const arma::mat VdualIon,const arma::cx_mat Ui,const arma::mat n);
arma::cx_mat solvePoisson(const operatorStruct O,const paramStruct P,arma::cx_mat n);
int getPsi(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi,
const arma::mat Vdual,std::shared_ptr<arma::cx_mat> returnPsi,std::shared_ptr<arma::mat> returnEpsilon);
double Prod(const arma::cx_mat A,const arma::cx_mat B);
arma::cx_mat diagOuter(paramStruct Pa,const arma::cx_mat A,const arma::cx_mat B);
arma::mat getVdual(const operatorStruct O,const arma::mat V);
arma::cx_mat USquaredInverse(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi);



string getELatex();
string getgradLatex();
string getPsiLatex();
string computeDensityLatex();
string computeUm1();


inline arma::mat getVdual(const operatorStruct O,const arma::mat V)
{

//! \brief function calculates the dual potential as used in energy functional formulation: J+(O(J(V)))
//! \param O...operator struct
//! \param V...potential matrix
//! \return matrix with dual potential

arma::mat returnMatrix;
arma::cx_mat cMatrix;
cMatrix=*O.Jd*(*O.O*(*O.J*V));
myFunctions::cassert(cMatrix.is_finite(),true,"getVdual: element in cMatrix matrix is infinite",__FILE__,__LINE__);
myFunctions::cassert(!cMatrix.has_nan(),true,"getVdual: element in cMatrix matrix is nan",__FILE__,__LINE__);
returnMatrix=real(cMatrix);

return returnMatrix;
}

inline arma::cx_mat diagOuter(paramStruct Pa, const arma::cx_mat A,const arma::cx_mat B)
{
//! \brief function calculates product of two matrices, then sums over rows: sum_row(A c(B))
//! \param A: matrix A
//! \param B: matrix B
//! definition of arma::sum(mat,dim) with dim=0 sum over rows and dim=1 over columns

verbosity(Pa,"diagOuter: compute the diagonal of outer product",3,__FILE__,__LINE__);

cx_mat outMat;
outMat=arma::sum(A%conj(B),1); // sum over row elements
//! should be Sx1 matrix
return outMat;

}

inline string computeUm1()
{
string description=" compute inverse of wavefunc integral $\\vec{U}(\\mathbf{k})$: \\newline ";
description+= "$\\vec{U} = W^{\\dagger} \\mathbf{O} W$ ";
return description;
}

/*compute the inverse of WOW*/
inline arma::cx_mat Uinvers(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi)
{
//! \brief latex description of Uinvers operator
    return inv(Wi.t()*(*O.O*Wi));
}

/*compute the inverse of WOW and takes square*/
inline arma::cx_mat USquaredInverse(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi)
{
//! \brief latex description of Uinvers operator
    return sqrtmat(inv(Wi.t()*(*O.O*Wi)));
}

inline string computeDensityLatex()
{
//! \brief latex description of density operator
string description=" compute density at every space point $\\vec{n}(\\mathbf{r})$: \\newline ";
description+= "$\\vec{n} = diag((\\mathbf{I}WU^{-1})(\\mathbf{I}W)^{\\dagger})$ ";
return description;
}


/*compute the density from wavefunction coefficients*/
inline arma::mat computeDensityFromWavefuncs(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi,const arma::cx_mat Ui)
{

    //! \brief compute electronic density using wavefunctions
    //! \param O: operators in use in one struct
    //! \param P: parameter in use in one struct
    //! \param Wi: wavefunction describing electronic state
    //! \return density as Sk x 1 dimensional matrix

    verbosity(P,"computeDensityFromWavefuncs: return the density: f*int(W*W)",2,__FILE__,__LINE__);
    verbosity(P,"computeDensityFromWavefuncs: use the diagouter routine to save memory",2,__FILE__,__LINE__);
    cx_mat a=(*O.I*Wi)*Ui;
    cx_mat b=(*O.I*Wi);
    verbosity(P,"computeDensityFromWavefuncs: diagOuter computes sum_row(a*bdagger)",2,__FILE__,__LINE__);
    cx_mat creturnM=diagOuter(P,a,b);
    verbosity(P,"computeDensityFromWavefuncs: we get a complex part of size"+std::to_string(abs(accu(imag(creturnM)))),3,__FILE__,__LINE__);
    verbosity(P,"computeDensityFromWavefuncs: we get a real part of size"+std::to_string(abs(accu(real(creturnM)))),3,__FILE__,__LINE__);
    mat returnM=P.f*real(creturnM);
    verbosity(P,"computeDensityFromWavefuncs: density return matrix has size "+std::to_string(returnM.n_rows)+" x "+std::to_string(returnM.n_cols),3,__FILE__,__LINE__);
    myFunctions::cassert(returnM.is_finite(),true,"element in density matrix is infinite",__FILE__,__LINE__);
    myFunctions::cassert(!returnM.has_nan(),true,"element in density matrix is nan",__FILE__,__LINE__);
    return returnM;
}

/*compute the potential from density*/
inline arma::cx_mat solvePoisson(const operatorStruct O,const paramStruct P,const arma::mat ni)
{
    //! \brief solve the poisson equation from electron density distribution (electron potential phi from density)
    //! \param O: operators in use in one struct
    //! \param P: parameter in use in one struct
    //! \param ni: real matrix of density Sx1
    //! \return phi: real matrix of potential phi Sx1

    return (*O.Li*(-4.*arma::datum::pi*(*O.O*(*O.J*ni))));

}

/*compute the hamiltonian*/
inline arma::cx_mat HW(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi,const arma::mat VdualIon,
const arma::cx_mat Ui,const arma::mat n)
{

    //! \brief solve the hamiltonian operator
    //! \param O: operators in use stored in struct
    //! \param P: parameter in use stored in struct
    //! \param Wi: complex matrix of N wavefunctions in N columns
    //! \param VdualIon: dual ion potential as real one column matrix
    //! \param Ui: inverse of U matrix, U=W+OW
    //! \param n: one column matrix holding the density in real space

    verbosity(P,"HW: Veff",2,__FILE__,__LINE__);
    arma::cx_mat Veff(VdualIon.n_rows,1,fill::zeros);
    #ifdef CALC_EION
    verbosity(P,"HW: add ionic potential to Veff",2,__FILE__,__LINE__);
    Veff.set_real(VdualIon);
    #endif
    verbosity(P,"HW: Veff n_rows: "+std::to_string(Veff.n_rows),2,__FILE__,__LINE__);
    verbosity(P,"HW: Veff n_cols: "+std::to_string(Veff.n_cols),2,__FILE__,__LINE__);
    #ifdef CALC_VEE
    verbosity(P,"HW: add electron potential part to Vdual",2,__FILE__,__LINE__);
    arma::cx_mat phi=solvePoisson(O,P,n);
    verbosity(P,"HW: phi n_rows: "+std::to_string(phi.n_rows),2,__FILE__,__LINE__);
    verbosity(P,"HW: phi n_cols: "+std::to_string(phi.n_cols),2,__FILE__,__LINE__);
    Veff+=*O.Jd*(*O.O*phi);
    verbosity(P,"HW: add exchange potential to Veffective",2,__FILE__,__LINE__);
    #endif
    #ifdef CALC_EXC
    verbosity(P,"H: calculate Exc potential",2,__FILE__,__LINE__);
    arma::mat Exc=excVWN(P,n);
    verbosity(P,"H: Exc n_rows: "+std::to_string(Exc.n_rows),2,__FILE__,__LINE__);
    verbosity(P,"H: Exc n_cols: "+std::to_string(Exc.n_cols),2,__FILE__,__LINE__);

    myFunctions::cassert(Exc.is_finite(),true,"HW: element in dExc matrix is infinite",__FILE__,__LINE__);
    myFunctions::cassert(!Exc.has_nan(),true,"HW: element in dExc matrix is nan",__FILE__,__LINE__);

    verbosity(P,"H: compute exchange correlation potential contribution",2,__FILE__,__LINE__);

    Veff+=*O.Jd*(*O.O*(*O.J*Exc));
    #endif
    verbosity(P,"H: n n_rows: "+std::to_string(n.n_rows),2,__FILE__,__LINE__);
    verbosity(P,"H: n n_cols: "+std::to_string(n.n_cols),2,__FILE__,__LINE__);
    #ifdef CALC_DEXC
    verbosity(P,"H: calculate dExc the derivative of Exc potential",2,__FILE__,__LINE__);
    arma::mat dExc=excpVWN(n);

    myFunctions::cassert(dExc.is_finite(),true,"HW: element in dExc matrix is infinite",__FILE__,__LINE__);
    myFunctions::cassert(!dExc.has_nan(),true,"HW: element in dExc matrix is nan",__FILE__,__LINE__);

    verbosity(P,"H: dExc n_rows: "+std::to_string(dExc.n_rows),2,__FILE__,__LINE__);
    verbosity(P,"H: dExc n_cols: "+std::to_string(dExc.n_cols),2,__FILE__,__LINE__);
    verbosity(P,"H: add exchange potential gradient contribution to Veffective",2,__FILE__,__LINE__);
    Veff+=arma::diagmat(dExc)*(*O.Jd*(*O.O*(*O.J*n)));
    #endif
    //arma::cx_mat Veff = Vdual+*O.Jd*(*O.O*phi)+*O.Jd*(*O.O*(*O.J*Exc)) + arma::diagmat(dExc)*(*O.Jd*(*O.O*(*O.J*n)));
    verbosity(P,"H: assemble and return H*W",2,__FILE__,__LINE__);
    cx_mat returnMatrix;
        verbosity(P,"H:  now compute the hamiltonian ",2,__FILE__,__LINE__);
    #ifdef CALC_KIN_ONLY
    returnMatrix= -0.5*(*O.L*Wi);
    #else
        #ifdef CALC_KIN
        verbosity(P,"H:  assemble kinetic and potential parts ",2,__FILE__,__LINE__);
        returnMatrix= -0.5*(*O.L*Wi)+*O.Id*(arma::diagmat(Veff)*(*O.I*Wi));
        #else

        returnMatrix= *O.Id*(arma::diagmat(Veff)*(*O.I*Wi));
        //returnMatrix= *O.Id*(Veff*(*O.I*Wi));
        #endif
    #endif

    myFunctions::cassert(returnMatrix.is_finite(),true,"HW: element in HW matrix is infinite",__FILE__,__LINE__);
    myFunctions::cassert(!returnMatrix.has_nan(),true,"HW: element in HW matrix is nan",__FILE__,__LINE__);

    return returnMatrix;

}

inline string getgradLatex()
{

//! \brief holds latex description o gradient operator

string description=" compute $\\nabla_{W} E$: ";
description+= " $\\nabla_{W} E = f*(H-O(W(U^{-1}(W^{T}H))))U^{-1}$ ";
return description;
}


inline arma::cx_mat getgrad(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi,const arma::mat VdualIon)
{

    //! \brief get the gradient of the energy along W
    //! \param O: operators in use in a struct
    //! \param P: parameter in use in a struct
    //! \param Wi: wavefunctions describing electronic state
    //! \param VdualIon: ion potential in real space
    //! \return gradient of energy functional along W
    //!
    //! compute $\nabla_{W} E$ as follows:
    //! $\nabla_{W} E = f*(H-O(W(U^{-1}(W^{T}H))))U^{-1}$
    //! compute $U^{-1}$, then hamiltonian $H(W)$, then assemble the gradient

    verbosity(P,"getgrad: compute Ui",2,__FILE__,__LINE__);
    arma::cx_mat Ui=Uinvers(O,P,Wi);
    verbosity(P,"getgrad: compute n",2,__FILE__,__LINE__);
    arma::mat n=computeDensityFromWavefuncs(O,P,Wi,Ui);
    verbosity(P,"getgrad: compute HW",2,__FILE__,__LINE__);
    arma::cx_mat HwM=HW(O,P,Wi,VdualIon,Ui,n);
    verbosity(P,"return the gradient",2,__FILE__,__LINE__);
    return P.f*(HwM-*O.O*(Wi*(Ui*(Wi.t()*HwM))))*Ui;
}


inline string getELatex()
{
    //! \brief return the latex description for the energy functional function


    string description=" $E=E_{kin}+E_{ee}+E_{exc}+E_{ion}$ \\newline"
    " $E_{kin}={f \\over 2} \\Sigma_{ii} (W^{T} (O (W U^{-1}))) $ with $ U^{-1} = (Wi^{T} Wi)^{-1}$ ... the kinetic energy of the electrons \\newline"
    " $E_{ee}={1 \\over 2} n_{e}(\\vec{x}) J^{\\dagger}(O(\\phi(\\vec{k})))$... the energy in the electronic potential field \\newline"
    " $E_{ion}= V_{ion}(\\vec{x}) n_{e}(\\vec{x}) ... the energy of the electrons in the ionic potential field$ \\newline";
    return description;
}

inline double getE(const operatorStruct O,const paramStruct P,const arma::cx_mat Wi,const arma::mat VdualIon)
{

    //! \brief compute the current energy of the electronic system
    //! \param O: operators in use in one struct
    //! \param P: parameter in use in one struct
    //! \param Wi: wavefunction describing electronic state
    //! \param VdualIon: ion potential in real space
    //! \return energy of current state as number
    //!
    //! $E=E_{kin}+E_{ee}+E_{exc}+E_{ion}$
    //! $E_{kin}={f \over 2} \Sigma_{ii} (W^{T} (O (W U^{-1}))) $ with $ U^{-1} = (Wi^{T}*Wi)^{-1}$ ... the kinetic energy of the electrons
    //! $E_{ee}={1 \over 2} n_{e}(\vec{x}) J^{\dagger}(O(\phi(\vec{k})))$... the energy in the electronic potential field
    //! $E_{ion}= V_{ion}(\vec{x}) n_{e}(\vec{x}) ... the energy of the electrons in the ionic potential field$


    verbosity(P,"getE: first compute the inverse of the wf",2,__FILE__,__LINE__);
    arma::cx_mat Ui=Uinvers(O,P,Wi);
    verbosity(P,"getE: then compute the density using the wf",2,__FILE__,__LINE__);
    arma::mat nm=computeDensityFromWavefuncs(O,P,Wi,Ui);
    //arma::cx_mat testM(*O.J*nm);
    //arma::cx_mat phi(*O.Li*(-4.*PI*(*O.O*(*O.J*nm))));
    //verbosity(P,"getE: build the Veffective operator",2,__FILE__,__LINE__);
    std::complex<double> cReturnEnergy(0,0);
    //cReturnEnergy=0;
    #ifdef CALC_KIN
        verbosity(P,"getE: first the kinetic part",2,__FILE__,__LINE__);
        cReturnEnergy=(-P.f*0.5*arma::trace(Wi.t()*(*O.L*(Wi*Ui))));
    #endif
    #ifdef CALC_VEE
    verbosity(P,"getE: next the electron potential self energy",2,__FILE__,__LINE__);
    verbosity(P,"getE: we solve the poisson equation",2,__FILE__,__LINE__);
    arma::cx_mat phi=solvePoisson(O,P,nm);
    arma::cx_mat tmpMat(phi);
    tmpMat=*O.O*phi;
    verbosity(P,"getE: compute the back transform to real space",2,__FILE__,__LINE__);
    tmpMat=*O.Jd*tmpMat;
    verbosity(P,"getE: now do dot product and transform to scalar",2,__FILE__,__LINE__);
    //cReturnEnergy+=(as_scalar(0.5*nm.diag().t()*tmpMat.diag()));
    cReturnEnergy+=(as_scalar(0.5*dot(nm.t(),tmpMat)));
    //returnEnergy=returnEnergy+real(arma::as_scalar(0.5*nm.t()*(*O.Jd*(*O.O*phi))));
    #endif
    #ifdef CALC_EXC
    verbosity(P,"getE: finally the exchange correlation energy",2,__FILE__,__LINE__);
    verbosity(P,"getE: calculate the exch potential",2,__FILE__,__LINE__);
    arma::mat Exc=excVWN(P,nm);
    cReturnEnergy+=(arma::as_scalar(vectorise(nm).t()*vectorise(*O.Jd*(*O.O*(*O.J*Exc)))));
    //returnEnergy=returnEnergy+real(arma::as_scalar(nm.t()*(*O.Jd*(*O.O*(*O.J*Exc)))));
    #endif
    #ifdef CALC_EION
    verbosity(P,"getE: then the electrons in the ion potential",2,__FILE__,__LINE__);
    //verbosity(P,"getE: Vdualinputvec: nrows"+std::to_string(Vdualinputvec.n_rows),2,__FILE__,__LINE__);
    //verbosity(P,"getE: Vdualinputvec: ncols"+std::to_string(Vdualinputvec.n_cols),2,__FILE__,__LINE__);
    verbosity(P,"getE: nm: nrows"+std::to_string(nm.n_rows),2,__FILE__,__LINE__);
    verbosity(P,"getE: nm: ncols"+std::to_string(nm.n_cols),2,__FILE__,__LINE__);
    //cReturnEnergy+=(arma::as_scalar(dot(vectorise(Vdualinputvec),nm.diag())));
    verbosity(P,"getE: Vdualinputvec: add (Vvec+)*nvec",2,__FILE__,__LINE__);
    verbosity(P,"getgrad: compute Vdual",2,__FILE__,__LINE__);
    cReturnEnergy+=(arma::as_scalar(dot(VdualIon.t(),nm)));
    #endif
    verbosity(P,"getE: now return the energy value",2,__FILE__,__LINE__);
    //arma::mat tmp(Wi.t()*Wi);
    double returnEnergy=real(cReturnEnergy);

    verbosity(P,"getE: we get a complex part of size"+std::to_string(abs(accu(imag(cReturnEnergy)))),3,__FILE__,__LINE__);
    verbosity(P,"getE: we get a real part of size"+std::to_string(abs(accu(real(cReturnEnergy)))),3,__FILE__,__LINE__);

    return returnEnergy;
}



inline string getPsiLatex()
{
    //! \brief return the latex description for the psi function

    string description=" $ \\Psi = Y D$ with $Y=W U^{-1/2}$ and $D: (DY)^{\\dagger} H (DY) = diag(\\vec{\\epsilon})$ \\newline";
    return description;
}


inline int getPsi(const operatorStruct O,const paramStruct P,const arma::cx_mat Win,const arma::mat Vdual,std::shared_ptr<arma::cx_mat> returnPsi
,std::shared_ptr<arma::mat> returnEpsilon)
{

    //! \brief calculate Psi, the final solution and real wavefunction
    //!
    //! \param O: operators in use in one struct
    //! \param P: parameter in use in one struct
    //! \param Win: wavefunction describing electronic state
    //! \param Vdual: ion potential in real space
    //! \param returnPsi: $\Psi$ returned to argument pointer to complex matrix
    //! \param returnEpsilon: energies of the various eigenfunctions of the Hamiltonian as real vector
    //! \return 1 for success
    //! with Y defined as $Y=W U^{-1/2}$
    //! $\Psi = Y D$ with $D: (DY)^{\dagger} H (DY) = diag(\vec{\epsilon})$


   // cout << arma::inv(arma::sqrt(Win.t()*(*O.O*Win)))*(arma::sqrt(Win.t()*(*O.O*Win))) << endl;
    //arma::cx_mat Ui=arma::inv(arma::sqrt(Win.t()*(*O.O*Win))); //Uinvers(O,P,Win);

    arma::cx_mat Ui=Uinvers(O,P,Win);
    arma::cx_mat Y=Win*USquaredInverse(O,P,Win); //*sqrt(Ui);
    arma::mat n=computeDensityFromWavefuncs(O,P,Win,Ui);
    arma::cx_mat HY=HW(O,P,Y,Vdual,Ui,n);
    verbosity(P,"getPsi: compute mu",2,__FILE__,__LINE__);
    arma::cx_mat mu=Y.t()*HY;
    //mu is hermitian -> mu_dag=mu
    cout << "mu " << endl;
    cout << mu-mu.t() << endl;

    arma::cx_vec epsilonVector(P.number_of_wavefunctions);

        for(int i=0;i<epsilonVector.n_elem;i++)
    {
        verbosity(P,"getPsi: after ini computation - eigenvalue "+std::to_string(i)+" is "+std::to_string(real(epsilonVector(i)))+" "+std::to_string(imag(epsilonVector(i))),2,__FILE__,__LINE__);
    }



    verbosity(P,"getPsi: compute eigenenergies and vectors",2,__FILE__,__LINE__);
    arma::eig_gen(epsilonVector,*returnPsi,mu);//solve for eigenvalues and eigenvectors

    #ifdef CHECK_NAN_FINITE

        myFunctions::cassert(returnPsi->is_finite(),ISCRITICAL,"getPsi: element in returnPsi matrix is infinite",__FILE__,__LINE__);
        myFunctions::cassert(!returnPsi->has_nan(),ISCRITICAL,"getPsi: element in returnPsi matrix is nan",__FILE__,__LINE__);
        myFunctions::cassert(epsilonVector.is_finite(),ISCRITICAL,"getPsi: element in epsilonVector matrix is infinite",__FILE__,__LINE__);
        myFunctions::cassert(!epsilonVector.has_nan(),ISCRITICAL,"getPsi: element in epsilonVector matrix is nan",__FILE__,__LINE__);
    #endif

    for(int i=0;i<epsilonVector.n_elem;i++)
    {
        verbosity(P,"getPsi: eigenvalue "+std::to_string(i)+" is "+std::to_string(real(epsilonVector(i)))+" "+std::to_string(imag(epsilonVector(i))),2,__FILE__,__LINE__);
    }

    verbosity(P,"getPsi: return values",2,__FILE__,__LINE__);
    returnEpsilon->col(0)=arma::real(epsilonVector);


    *returnPsi=Y*(*returnPsi);

    cout << Y.t()*Y << endl;
    cout << "should be the identity: " << endl;
    cout << returnPsi->t()*(*O.O*(*returnPsi)) << endl;
    cout << "should be diagonal: " << endl;
    cout << returnPsi->t()*HW(O,P,*returnPsi,Vdual,Ui,n) << endl;
    cout << "should match the following vector entries: " << endl;
    cout << epsilonVector << endl;

    return 1;
}

inline double Prod(const arma::cx_mat A,const arma::cx_mat B)
{
// \brief memory saving method of building product of A and B
    return real(trace(A.t()*B));
}

#endif // DFTFUNCTIONS_H_INCLUDED
