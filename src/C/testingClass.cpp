#include "testingClass.h"
#include "slice.h"
#include "smooth.h"

using namespace myFunctions;
using namespace arma;

testingClass::testingClass(latexComment *latX,paramStruct P)
{
    _myLatexPtr=latX;
    _myP=P;
}

testingClass::~testingClass()
{
    //dtor
}


bool testingClass::schroedingerTest(const operatorStruct Op,const paramStruct Pa,arma::Col<double> S,Mat<double> R,mat r,
mat mat_center_of_cell,const arma::cx_mat Sf,arma::mat X,arma::Mat<double> G2)
{
  //! \brief test for the solution of the Schroedinger equation and a parabolic model potential
    //!
    //! \param O an operator struct
    //! \param P a parameter struct
    //! \param Sgrid - matrix with 3 entries, with number of points along each direction in space
    //! \param r coordinates of cell Sx3 matrix
    //! \param mat_center_of_cell: center of cell coordinates, Sx3 matrix
    //! \param X: coordinates of atom in system

    //! \return true if test passed, else false

        _myLatexPtr->subSection("Solve the Schroedinger equation using an oscillator potential");
        _myLatexPtr->subsubSection("define an oscillator potential");
        _myLatexPtr->newLine(" $V(\\mathbf{r})={1 \\over 2} \\omega^{2} d\\mathbf{r}^{2}$ ");
        verbosity(_myP,"schroedinger: solve the schroedinger equation",2,__FILE__,__LINE__);
        verbosity(_myP,"schroedinger: define the oscillator potential",2,__FILE__,__LINE__);
        #define SCHROEDINGER_PREPROCESSING
        #ifdef SCHROEDINGER_PREPROCESSING
        double omega=2.;
        mat dr(arma::sqrt(arma::sum((r-mat_center_of_cell)%(r-mat_center_of_cell),1)));
        arma::mat Vosci(dr.n_rows,1);
        Vosci=0.5*pow(omega,2)*(dr%dr);

        gnuPlotPlotting gpL;
        string imgName=gpL.plotAMatrixSlice(_myP,"schroedinger_Vosci",Vosci,Pa.S,0);
        _myLatexPtr->newLine("A parabolic potential is assumed with the following shape:");
        _myLatexPtr->insertImage(imgName,"parabolic ion potential for which Schroedinger equation is solved.",Pa.caseName);

        verbosity(_myP,"schroedinger: the oscillator potential has dimension "+std::to_string(Vosci.n_rows),2,__FILE__,__LINE__);
        verbosity(_myP,"schroedinger: initialize wavefunction W",2,__FILE__,__LINE__);
        _myLatexPtr->subsubSection("Preparation of the wavefunction");
        //arma::cx_mat Wm(prodS,Pa.number_of_wavefunctions);
        _myLatexPtr->newLine(" Initialize wavefunctions randomĺy : ");
        _myLatexPtr->newLine(" $W_{m}:$ complex matrix with $ \\Sigma S_{k} \\times N_{s} $ elements in r space ");

        //Wm=arma::randu<arma::cx_mat>(prodS,Pa.number_of_wavefunctions);
        //std::shared_ptr<arma::cx_mat> W(&Wm);
        std::shared_ptr<arma::cx_mat> W(new arma::cx_mat(arma::randn<arma::cx_mat>(Pa.prodS,Pa.number_of_wavefunctions)));

        #endif
        #define SCHROEDINGER_DO_SD
        #ifdef SCHROEDINGER_DO_SD
        verbosity(_myP,"schroedinger: create optimizer sd",2,__FILE__,__LINE__);
        arma::mat Vdual=getVdual(Op,Vosci);
        sdOptimizer sd(Op,Pa,Vdual,_myLatexPtr);

        sd.myLatexClass->commentMyFunction();
        verbosity(_myP,"schroedinger: prepare sd optimizer",2,__FILE__,__LINE__);
        sd.setup(W);
        verbosity(_myP,"schroedinger: compute density and plot",2,__FILE__,__LINE__);

        //check if W is normalized W=Y=W*U^⁻1/2
        verbosity(_myP,"schroedinger: check if W is normalized W=Y=W*U^⁻1/2",2,__FILE__,__LINE__);

        cx_mat Ui=Uinvers(Op,Pa,*W.get());
        mat density(computeDensityFromWavefuncs(Op,Pa,*W.get(),Ui));
        gpL.plotMatrix3Slices(_myP,"schroedinger_initialDensity",density,Pa.S);

        _myLatexPtr->newLine("The following image shows the initial density computed using the initial state of random wavefunction.");

        _myLatexPtr->insertImage("schroedinger_initialDensity_xy","density after initialization of the wavefunctions",Pa.caseName);

        sd.myLatexClass->commentMyFunction();
        verbosity(_myP,"====================================",2,__FILE__,__LINE__);
        verbosity(_myP,"schroedinger: optimize using sd",2,__FILE__,__LINE__);
        verbosity(_myP,"====================================",2,__FILE__,__LINE__);

        sd.optimize(50,W);

        pccgOptimizer pcg(Op,Pa,Vdual,G2,_myLatexPtr);
        pcg.orthogonalizeWfunc(W);
        double finalEnergy=pcg.optimize(250,W);

        sd.myLatexClass->commentMyFunction();
        verbosity(_myP,"now get the wavefunction",2,__FILE__,__LINE__);
        // Psi as a pointer to the wavefunction
        std::shared_ptr<arma::cx_mat> Psi(new arma::cx_mat(*W));
        // Epsilon points to the eigenvalues (real vector)
        //std::shared_ptr<arma::cx_vec> Epsilon(new arma::cx_vec(W->n_rows));
        std::shared_ptr<arma::mat> Epsilon(new arma::mat(Pa.number_of_wavefunctions,1));

        #endif // SCHROEDINGER_DO_SD
        #define SCHROEDINGER_CALC_PSI
        #ifdef SCHROEDINGER_CALC_PSI
        if(getPsi(Op,Pa,*W,Vosci,Psi,Epsilon))
        {
            cout << "succesfully extracted eigenvalues and eigenvectors from the final solution W!" << std::endl;
             cout << "schroedinger: returnEpsilon " << Epsilon->col(0) << endl;

        }
        _myLatexPtr->newLine("The electron states have the following energies:");

        #endif
        #define SCHROEDINGER_PLOT_RESULT
        #ifdef SCHROEDINGER_PLOT_RESULT
        arma::mat dat(Psi->n_rows,1);
        string name;
        for(unsigned int st=0; st<Psi->n_cols; st++)
        {
             cout << "schroedinger: returnEpsilon " << Epsilon->row(0) << endl;

            //arma::cx_mat AEpsilon = cx_mat(Epsilon.get()->zeros(),Epsilon.get()->zeros());
            std::cout << "=== State " << st << ", has energy = " << Epsilon->row(st) << " === " << std::endl;
            _myLatexPtr->newLine("=== state " + std::to_string(st) + " has energy " + std::to_string(as_scalar(Epsilon->row(st)))+ "hartree ===");

            verbosity(_myP,"now we get the square of the wavefunction for output",2,__FILE__,__LINE__);
            dat=real(pow(*Op.I*Psi->col(st),2));
            verbosity(_myP,"print three slices of result to as ppm / gnuplot",2,__FILE__,__LINE__);
            #ifdef PLOT_RESULT
            for(int k=0; k<3; k++)
            {
                string plane;
                if(k==0) {plane="yz";}
                else if(k==1) {plane="xz";}
                else if(k==2) {plane="xy";}
                verbosity(_myP,"now take a slice of data at plane "+plane,2,__FILE__,__LINE__);
                arma::mat sl;
                verbosity(_myP,"take a slice "+plane,2,__FILE__,__LINE__);
                sl=myFunctions::slice(Pa,dat,Pa.S,(int) Pa.S(k)/2.,k);
                #ifdef PLOT_PPM
                verbosity(_myP,"get name for ppm file "+plane,2,__FILE__,__LINE__);
                name="psi"+std::to_string(st)+"d_m_"+plane; //std::to_string(k)
                verbosity(_myP,"write ppt file "+name,2,__FILE__,__LINE__);
                ppm(name,sl*0.3,sl,sl,&latX);
                #endif
                #ifdef PLOT_GNUPLOT
                if(st==0) _myLatexPtr->newLine("The physical eingenstates of the hamiltonian Psi are shown in the next figure:");
                if(k==0)
                {
                    name="schroedinger_psi_"+std::to_string(st);
                    string imageName=gpL.plotAMatrixSlice(Pa,name,dat,Pa.S,k);
                    string caption="result from solution of schroedinger equation: slice through wave function of electron " + std::to_string(st) + " plane "+plane;
                    _myLatexPtr->insertImage(imageName,caption,Pa.caseName);
                }
                #endif
            }
            #endif
        }
        #endif

        Ui=Uinvers(Op,Pa,*W.get());
        density=computeDensityFromWavefuncs(Op,Pa,*W.get(),Ui);
        imgName=gpL.plotAMatrixSlice(Pa,"schroedinger_density_final",density,Pa.S,0);
        _myLatexPtr->newLine("Finally a density distribution as shown in the next image is computed.");
        _myLatexPtr->insertImage(imgName,"final result of Schroedinger equation: density",Pa.caseName);

        bool testPassed=false;
        if(abs(finalEnergy-18)<0.01)
            {
            testPassed=true;
            _myLatexPtr->newLine(" Passed the Schroedinger equation test ");
            }
            else
            {
            _myLatexPtr->newLine(" Failed the Schroedinger equation test ");
            }
        return testPassed;


}

bool testingClass::multicolumnTest(const operatorStruct O,const paramStruct P,arma::Col<double> S,
Mat<double> R,mat r,mat mat_center_of_cell,const arma::cx_mat Sf,arma::mat X)
{
  //! \brief check for multicolumn consistency of operators
    //!
    //! \param O an operator struct
    //! \param P a parameter struct
    //! \param Sgrid - matrix with 3 entries, with number of points along each direction in space
    //! \param r coordinates of cell Sx3 matrix
    //! \param mat_center_of_cell: center of cell coordinates, Sx3 matrix
    //! \param X: coordinates of atom in system

    //! \return true if test passed, else false

    bool testPassed=false;

    cx_mat in(prod(S),3);
    in.set_real(randn<mat>(prod(S),3));
    in.set_imag(randn<mat>(prod(S),3));
    cx_mat outM=*O.I*in;
    cx_mat outS(outM);
    cx_mat temp=*O.I*in.col(0);//.submat(0,in.n_rows,0,0);
    outS.col(0)=temp.col(0);
    outS.col(1)=*O.I*in.col(1);
    outS.col(2)=*O.I*in.col(2);
    std::complex<double> diff=(abs(outM-outS)).max();

    const std::complex<double> bound(1.0e-5,1.0e-05);
    if(std::real(diff*std::conj(diff))<std::real(bound*std::conj(bound))){testPassed=true;}
    return testPassed;

}

bool testingClass::checkHermitianConsistency(const operatorStruct O,const paramStruct P,arma::Col<double> S,
Mat<double> R,mat r,mat mat_center_of_cell,const arma::cx_mat Sf,arma::mat X)
{
  //! \brief test hermitian operators
    //!
    //! \param O an operator struct
    //! \param P a parameter struct
    //! \param Sgrid - matrix with 3 entries, with number of points along each direction in space
    //! \param r coordinates of cell Sx3 matrix
    //! \param mat_center_of_cell: center of cell coordinates, Sx3 matrix
    //! \param X: coordinates of atom in system

    //! \return true if test passed, else false

    bool testPassed=false;

    cx_mat a(prod(S),1);
    cx_mat b(prod(S),1);
    a.set_real(randn<mat>(prod(S),1));
    b.set_real(randn<mat>(prod(S),1));
    a.set_imag(randn<mat>(prod(S),1));
    b.set_imag(randn<mat>(prod(S),1));
    std::complex<double> diff=as_scalar(sum(conj(a.t()*(*O.I*b))-b.t()*(*O.Id*a)));
    diff+=as_scalar(sum(conj(a.t()*(*O.J*b))-b.t()*(*O.Jd*a)));

    const std::complex<double> bound(1.0e-5,1.0e-05);
    if(std::real(diff*std::conj(diff))<std::real(bound*std::conj(bound))){testPassed=true;}
    return testPassed;


}


bool testingClass::ewaldTest(const operatorStruct O,const paramStruct P,arma::Col<double> S,
Mat<double> R,mat r,mat mat_center_of_cell,const arma::cx_mat SfM,arma::mat X)
{
  //! \brief test ewald sums of a simple physical system
    //!
    //! \param O an operator struct
    //! \param P a parameter struct
    //! \param Sgrid - matrix with 3 entries, with number of points along each direction in space
    //! \param r coordinates of cell Sx3 matrix
    //! \param mat_center_of_cell: center of cell coordinates, Sx3 matrix
    //! \param X: coordinates of atom in system

    //! \return true if test passed, else false

_myLatexPtr->subSection("The ewald summation test");
_myLatexPtr->newLine(" We solve the following equation and compare with an analytical result: ");
_myLatexPtr->newLine(" $\\nabla^{2} \\phi = n $ }");
_myLatexPtr->newLine(" In operator language this writes as: ");
_myLatexPtr->newLine(" $I(L(-4 \\pi ( O (J(n)))))$ ");
_myLatexPtr->newLine(" J: $n(r) \\mapsto n(k)$ ");
_myLatexPtr->newLine(" O: integral over product of basis vectors from $L^{-1},J(n)$ ");
_myLatexPtr->newLine(" L: the inverse Laplace transform in k - space - $ n(k) \\mapto {n(k) \\over k^{2}} $ ");
_myLatexPtr->newLine(" I: ${n(k) \\over k^{2}} \\mapsto n(k)$  ");


_myLatexPtr->newLine(" compute distance to center as: $dr = \\sqrt{ \\Sigma_{i} (r-\\langle r \\rangle)^{2} $ }");


Mat<double> dr(arma::sqrt(arma::sum((r-mat_center_of_cell)%(r-mat_center_of_cell),1)));


const double sigma1=0.25;



_myLatexPtr->newLine(" charge distribution 1 as: $g_{1}=Z{\\exp(-{dr^{2} \\over 2 \sigma_{1})^{2}}) \\over 2 \\pi \\sigma_{1}^{3}}$");

verbosity(_myP,"compute g1",2,__FILE__,__LINE__);

arma::cx_mat g1(dr.n_rows,1,arma::fill::zeros);
g1.set_real(P.Z*arma::exp(-(dr%dr)/(2.*pow(sigma1,2)))/sqrt(pow(2.*arma::datum::pi*sigma1*sigma1,3)));

verbosity(_myP,"sum dr = "+std::to_string(arma::as_scalar(arma::sum(arma::sum(dr)))),2,__FILE__,__LINE__);
verbosity(_myP,"sum g1 = "+std::to_string(arma::as_scalar(real(arma::sum(arma::sum(g1))))),2,__FILE__,__LINE__);

cassert(!dr.has_inf(),ISCRITICAL,"dr has inf",__FILE__,__LINE__);
cassert(!dr.has_nan(),ISCRITICAL,"dr has nan",__FILE__,__LINE__);

cassert(!g1.has_inf(),ISCRITICAL,"g1 has inf",__FILE__,__LINE__);
cassert(!g1.has_nan(),ISCRITICAL,"g1 has nan",__FILE__,__LINE__);


//arma::cx_mat n(g1.n_elem,1);
//n.fill(0);
verbosity(_myP,"now extract density by applying the charge distribution to all the fields using the structure function Sf",2,__FILE__,__LINE__);

_myLatexPtr->newLine(" charge density n is $n=g_{2}-g_{1}$");

arma::cx_mat nc(arma::size(dr),fill::zeros);

cassert(!SfM.has_inf(),ISCRITICAL,"Sf has inf",__FILE__,__LINE__);
cassert(!SfM.has_nan(),ISCRITICAL,"Sf has nan",__FILE__,__LINE__);
cassert(SfM.is_finite(),ISCRITICAL,"Sf has nan",__FILE__,__LINE__);

nc.set_real(real(*O.I*((*O.J*g1)%SfM)));

cassert(!nc.has_inf(),ISCRITICAL,"nc has inf",__FILE__,__LINE__);
cassert(!nc.has_nan(),ISCRITICAL,"nc has nan",__FILE__,__LINE__);
cassert(nc.is_finite(),ISCRITICAL,"nc has nan",__FILE__,__LINE__);

//gpL.plotMatrix3Slices("ewaldrealSf",real(SfM),S);
//gpL.plotMatrix3Slices("ewaldimagSf",imag(SfM),S);


arma::mat nreal=real(nc);

//gpL.plotMatrix3Slices("nreal",nreal,S);

verbosity(_myP,"ewald: compute density distribution using structure function Sf",2,__FILE__,__LINE__);

verbosity(_myP,"max nc = "+std::to_string(arma::as_scalar(real(arma::max(arma::max(nc))))),2,__FILE__,__LINE__);
verbosity(_myP,"min nc = "+std::to_string(arma::as_scalar(real(arma::min(arma::min(nc))))),2,__FILE__,__LINE__);

verbosity(_myP,"max nreal = "+std::to_string(arma::as_scalar(arma::max(arma::max(nreal)))),2,__FILE__,__LINE__);
verbosity(_myP,"min nreal = "+std::to_string(arma::as_scalar(arma::min(arma::min(nreal)))),2,__FILE__,__LINE__);

////for(int dir=0;dir<3;dir++)
////{
////    verbosity(_myP,"ewald: nreal nrow "+std::to_string(nreal.n_rows)+" ncols "+std::to_string(nreal.n_cols),2,__FILE__,__LINE__);
////
////    TmpMat=myFunctions::slice(nreal,S,S(dir)/2.,dir);
////
////    arma::mat SmoothedMat;
////    SmoothedMat = smooth(TmpMat);
////
////    cassert(!SmoothedMat.has_inf(),isCritical,"SmoothedMat has inf",__FILE__,__LINE__);
////    cassert(!SmoothedMat.has_nan(),isCritical,"SmoothedMat has nan",__FILE__,__LINE__);
////
////    string fileName;
////    fileName = "testOutputSmoothedMatrixEwald_" + std::to_string(dir) +".2dMat";
////    ofstream of(fileName);
////    of << TmpMat << endl;
////    of.close();
////
////    //TmpMat=sliceIt(nreal,S,S(dir)/2.,dir);
////    verbosity(_myP,"ewald: TmpMat nrow "+std::to_string(TmpMat.n_rows)+" ncols "+std::to_string(TmpMat.n_cols),2,__FILE__,__LINE__);
////    verbosity(_myP,"ewald: SmoothedMat nrow "+std::to_string(SmoothedMat.n_rows)+" ncols "+std::to_string(SmoothedMat.n_cols),2,__FILE__,__LINE__);
////
////    //#ifdef PLOT_WITH_CV
////    //    openCVPlotting* pL(new openCVPlotting);
////    //    pL->plotMatrix("g1",SmoothedMat);//arma::real(g1));
////    //    delete pL;
////    //#endif
////
////    cassert(!TmpMat.has_inf(),isCritical,"TmpMat has inf",__FILE__,__LINE__);
////    cassert(!TmpMat.has_nan(),isCritical,"TmpMat has nan",__FILE__,__LINE__);
////    cassert(TmpMat.is_finite(),isCritical,"TmpMat is infinite",__FILE__,__LINE__);
////
////    //verbosity(_myP,"ewald: get name for ppm file ",2,__FILE__,__LINE__);
////    //string name="n_ewald_"+std::to_string(dir)+"d_m"; //std::to_string(k)
////    //verbosity(_myP,"write ppt file "+name,2,__FILE__,__LINE__);
////    //
////    //ppm(name,TmpMat*0.3,TmpMat,TmpMat,_myLatexPtr);
////}
//%# Check norms and integral (should be near 1 and 0, respectively)
cout << "Normalization check on g1:" << sum(g1)*det(R)/prod(S) << endl;
cout << "Total charge check:" << sum(real(nc))*det(R)/prod(S) << endl;

//_myLatexPtr->newLine(" check for normalization 1 $ \\Sigma_{i} g_{i} {V \\over N} = \\Sigma_{i} g_{i} dV$ = "+std::to_string(sum(g1,0)*det(R)/prod(S)));
//_myLatexPtr->newLine(" with $V$ the volume of the computations mesh, N the number of grid points and dV the volume of one grid cell");

verbosity(_myP,"ewald: nc is a real stored in a complex matrix",2,__FILE__,__LINE__);
verbosity(_myP,"ewald: now compute potential phi",2,__FILE__,__LINE__);

//arma::cx_mat phi=*O.I*(*O.Li*(-4.*arma::datum::pi*(*O.O*(*O.J*nc))));

arma::cx_mat phik = solvePoisson(O,P,nreal);
arma::cx_mat phi = *O.I*phik;

double Unum=0.5*real(arma::accu(((*O.J*phi).t()*(*O.O*(*O.J*nc)))));

arma::mat Xtranspose=X.t();
verbosity(_myP,"ewald: we have "+std::to_string(Xtranspose.n_rows)+" atoms in our system!",2,__FILE__,__LINE__);

double Uself=P.Z*P.Z/(2.*sqrt(arma::datum::pi))*(1./sigma1)*Xtranspose.n_rows;


//((1./sigma1+1./sigma2)/2.-std::sqrt(2.)/sqrt(pow(sigma1,2)+pow(sigma2,2)))/sqrt(arma::datum::pi);
//fprintf("Numeric, analytic Coulomb energy: %20.16f,%20.16f\n",Unum,Uanal);
std::cout << "numeric coulomb energy: " << Unum << " analytic coulomb energy " << Uself << endl;

_myLatexPtr->newLine(" The numeric coulomb energy is "+std::to_string(Unum));
_myLatexPtr->newLine(" The analytically computed energy is "+std::to_string(Uself));

bool testPassed=false;
if(abs(Uself-Unum-0.33)<0.01)
    {
    testPassed=true;
    _myLatexPtr->newLine(" Passed the poisson test ");
    }
    else
    {
    _myLatexPtr->newLine(" Failed the poisson test ");
    }
return testPassed;
}


bool testingClass::poissonEquationTest(const operatorStruct O,const paramStruct P,arma::Col<double> S,Mat<double> R,mat r,mat mat_center_of_cell)
{

    //! \brief solve the poisson equation and compare with analytical result
    //!
    //! \param O an operator struct
    //! \param P a parameter struct
    //! \param Sgrid - matrix with 3 entries, with number of points along each direction in space
    //! \param r coordinates of cell Sx3 matrix
    //! \param mat_center_of_cell: center of cell coordinates, Sx3 matrix

    //! \return true if test passed, else false

_myLatexPtr->subSection("The poisson equation test");
_myLatexPtr->newLine(" We solve the following equation and compare with an analytical result: ");
_myLatexPtr->newLine(" $\\nabla^{2} \\phi = n $ }");
_myLatexPtr->newLine(" In operator language this writes as: ");
_myLatexPtr->newLine(" $I(L(-4 \\pi ( O (J(n)))))$ ");
_myLatexPtr->newLine(" J: $n(r) \\mapsto n(k)$ ");
_myLatexPtr->newLine(" O: integral over product of basis vectors from $L^{-1},J(n)$ ");
_myLatexPtr->newLine(" L: the inverse Laplace transform in k - space - $ n(k) \\mapto {n(k) \\over k^{2}} $ ");
_myLatexPtr->newLine(" I: ${n(k) \\over k^{2}} \\mapsto n(k)$  ");


_myLatexPtr->newLine(" compute distance to center as: $dr = \\sqrt{ \\Sigma_{i} (r-\\langle r \\rangle)^{2} $ }");


Mat<double> dr(arma::sqrt(arma::sum((r-mat_center_of_cell)%(r-mat_center_of_cell),1)));

cout << "sizedr" << sum(dr) << endl;

const double sigma1=0.75;
const double sigma2=0.5;

verbosity(_myP,"compute g1,g2",2,__FILE__,__LINE__);


_myLatexPtr->newLine(" charge distribution 1 as: $g_{1}={\\exp(-{dr^{2} \\over 2 \sigma_{1})^{2}}) \\over 2 \\pi \\sigma_{1}^{3}}$");

//Col<double> g1=arma::exp(-(dr%dr)/(2.*pow(sigma1,2)))/sqrt(pow(2.*arma::datum::pi*sigma1*sigma1,3));
mat g1=arma::exp(-(dr%dr)/(2.*pow(sigma1,2)))/sqrt(pow(2.*arma::datum::pi*sigma1*sigma1,3));

_myLatexPtr->newLine(" charge distribution 2 as: $g_{2}={\\exp(-{dr^{2} \\over 2 \sigma_{2})^{2}}) \\over 2 \\pi \\sigma_{2}^{3}}$");

//Col<double> g2=arma::exp(-(dr%dr)/(2.*pow(sigma2,2)))/sqrt(pow(2.*arma::datum::pi*sigma2*sigma2,3));
mat g2=arma::exp(-(dr%dr)/(2.*pow(sigma2,2)))/sqrt(pow(2.*arma::datum::pi*sigma2*sigma2,3));

verbosity(_myP,"compute n",2,__FILE__,__LINE__);

//arma::cx_mat n(g1.n_elem,1);
//n.fill(0);

_myLatexPtr->newLine(" charge density n is $n=g_{2}-g_{1}$");

//for(int i=0;i<n.n_cols;i++){n.col(i)=cx_vec(g2-g1,arma::zeros<arma::vec>(g2.n_elem));} //arma::fill::zeros);}

arma::mat nr(g1.n_elem,1,fill::zeros);
nr=g2-g1;
arma::cx_mat phik=solvePoisson(O,P,nr);
arma::cx_mat phi = *O.I*phik;

//%# Check norms and integral (should be near 1 and 0, respectively)
cout << "Normalization check on g1:" << sum(g1)*det(R)/prod(S) << endl;
cout << "Normalization check on g2:" << sum(g2)*det(R)/prod(S) << endl;
cout << "Total charge check:" << sum(nr)*det(R)/prod(S) << endl;

_myLatexPtr->newLine(" check for normalization 1 $ \\Sigma_{i} g_{i} {V \\over N} = \\Sigma_{i} g_{i} dV$ = "+std::to_string(sum(g1.col(0))*det(R)/prod(S)));
_myLatexPtr->newLine(" check for normalization 2 $ \\Sigma_{i} g_{i} {V \\over N} = \\Sigma_{i} g_{i} dV$ = "+std::to_string(sum(g2.col(0))*det(R)/prod(S)));
_myLatexPtr->newLine(" with $V$ the volume of the computations mesh, N the number of grid points and dV the volume of one grid cell");

//arma::cx_mat phi=*O.I*(*O.Li*(-4.*arma::datum::pi*(*O.O*(*O.J*n))));


ofstream outp;
outp.open("poisson_phi.dat");
outp << phi << endl;
outp.close();


//cout << arma::accu(phi)<< endl;
//arma::cx_mat intermediateM(-4.*arma::datum::pi*(*O.O*(*O.J*n)));
//cout << arma::accu(intermediateM.rows(0,12)) << endl;
//exit(-1);
/*arma::cx_mat toaccu((J*phi).t()%(Oo*(J*n)));
std::complex<double> cxvart=arma::accu(toaccu);
double num=real(cxvart);*/
double Unum=0.5*real(arma::accu((*O.J*phi)%(*O.O*(*O.J*nr))));
double Uanal=((1./sigma1+1./sigma2)/2.-std::sqrt(2.)/sqrt(pow(sigma1,2)+pow(sigma2,2)))/sqrt(arma::datum::pi);
//fprintf("Numeric, analytic Coulomb energy: %20.16f,%20.16f\n",Unum,Uanal);
std::cout << "numeric coulomb energy: " << Unum << " analytic coulomb energy " << Uanal << endl;

_myLatexPtr->newLine(" The numeric coulomb energy is "+std::to_string(Unum));
_myLatexPtr->newLine(" The analytical coulomb energy is "+std::to_string(Uanal));

bool testPassed=false;
if(abs(Unum-Uanal)<0.01)
    {
    testPassed=true;
    _myLatexPtr->newLine(" Passed the poisson test ");
    }
    else
    {
    _myLatexPtr->newLine(" Failed the poisson test ");
    }
return testPassed;
}

bool testingClass::IJtest(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid)
{
    //! \brief do a J I inversion transformation test
    //!
    //! \param O an operator struct
    //! \param P a parameter struct
    //! \return S, the S column vector containing the grid dimension

bool returnval=false;
arma::cx_mat dW=arma::randu<arma::cx_mat>(arma::prod(Sgrid),1);
//cout << dW << endl;
arma::cx_mat result(*O.J*(*O.I*dW)/dW);
//cout << real(sum(*O.J*(*O.I*dW)-dW));
double realdifference=as_scalar(real(sum(*O.J*(*O.I*dW)-dW)));
double imagdifference=as_scalar(imag(sum(*O.J*(*O.I*dW)-dW)));

if((realdifference+imagdifference)<0.01){returnval=true;}
return returnval;
}

bool testingClass::IJtestMultiColumn(const operatorStruct O,const paramStruct P,arma::Col<double> Sgrid)
{
    //! \brief do a J I inversion transformation test with a multicolumn matrix input
    //!
    //! \param O an operator struct
    //! \param P a parameter struct
    //! \return S, the S column vector containing the grid dimension

bool returnval=false;
arma::cx_mat dW=arma::randu<arma::cx_mat>(arma::prod(Sgrid),5);
//cout << dW << endl;
arma::cx_mat result(*O.J*(*O.I*dW)/dW);
//cout << real(sum(*O.J*(*O.I*dW)-dW));
arma::vec rowSums=vectorise(real(sum(*O.J*(*O.I*dW)-dW)));
double realdifference=as_scalar(sum(rowSums,0));
rowSums=vectorise(imag(sum(*O.J*(*O.I*dW)-dW)));
double imagdifference=as_scalar(sum(rowSums,0));

if((realdifference+imagdifference)<0.01){returnval=true;}
return returnval;
}
