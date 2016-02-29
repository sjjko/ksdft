#include "atomicSystem.h"

using namespace std;
using namespace arma;


int atomicSystem::setupPotential()
{
 //! \brief set up the computational potential - output is the dual potential -
 //! Vdual lives in real space - enters the electron energy as follows: $E_{ion}=Vdual^{\dagger} n$
//! \return Vdual: matrix holds dual potential
//! uses source code definition which has to reside in the case dictionary
//! uses dynamic linking and compilation at run time

verbosity(_Pa,"atomicSystem::setupPotential get potential using ionicpotentialClass",2,__FILE__,__LINE__);
ionicPotentialClass pot(_Op,_Pa,_dr,_Sf,_X,_G2,_Pa.caseName);
pot.computePotential();
pot.checkPotential();
this->_Vdual=pot.getPotential();
verbosity(_Pa,"atomicSystem::setupPotential Vdual read",2,__FILE__,__LINE__);

//string pathToCaseDir=" ./case/"+_Pa.caseName+"/";
//string includeDirectories="-I /home/kos/Code/kos/minidft/ksdft++/case/"+_Pa.caseName+"/ ";
//includeDirectories+=" -I /home/kos/Code/kos/minidft/ksdft++/src/H/ ";
//includeDirectories+=" -I /home/kos/Code/kos/minidft/ksdft++/src/C/ ";
//includeDirectories+="-I /home/kos/Code/kos/minidft/ksdft++/texcaller/c/ ";
//string sourceFile="./src/C/dualPotential.cpp";
//verbosity(_Pa,"the main potential source file resides with the other headers"+sourceFile,0,__FILE__,__LINE__);
//verbosity(_Pa,"check if we can find the potential source file "+sourceFile,2,__FILE__,__LINE__);
//ifstream f(sourceFile.c_str());
//myFunctions::cassert(f.good(),ISCRITICAL,"could not find dualPotential.cpp file in above directory - Exit!",__FILE__,__LINE__);
//f.close();
//verbosity(_Pa,"compile the ionic potential source file "+sourceFile,0,__FILE__,__LINE__);
//string libraryFile="/home/kos/Code/kos/minidft/ksdft++/case/"+_Pa.caseName +"/dualPotential.so ";
//verbosity(_Pa,"remove existing library files "+sourceFile,0,__FILE__,__LINE__);
//string systemCallString="rm "+libraryFile;
////system(systemCallString.c_str());
//string externalLibraries=" /home/kos/Code/kos/minidft/ksdft++/texcaller/c/libtexcaller.so ";
////externalLibraries+=" ./obj/Debug/src/C/cJ.o ";
////externalLibraries+=" ./obj/Debug/src/C/cI.o ";
////externalLibraries+=" ./obj/Debug/src/C/cJdag.o ";
////externalLibraries+=" ./obj/Debug/src/C/cIdag.o ";
////externalLibraries+=" ./obj/Debug/src/C/Linv.o ";
////externalLibraries+=" ./obj/Debug/src/C/L.o ";
////externalLibraries+=" ./obj/Debug/src/C/O.o ";
//systemCallString="c++ -std=c++11 "+sourceFile+" " + includeDirectories + " -o ./dualPotential.so "+externalLibraries+" -shared -fPIC"; //-larmadillo -ldl
//cout << "We make the following system call: " << systemCallString << endl;
//system(systemCallString.c_str());
//verbosity(_Pa,"link the dualPotential library "+libraryFile,0,__FILE__,__LINE__);
//cout << libraryFile.c_str() << " open libraryFile.c_str() " << endl;
////void* potentialHandle = dlopen(libraryFile.c_str(), RTLD_LAZY); //RTLD_NOW);
//void* potentialHandle = dlopen("/home/kos/Code/kos/minidft/ksdft++/dualPotential.so",RTLD_NOW);
//cout << " finished dlopen" << endl;
////cout << dlerror() << " the dlerror " << endl;
////myFunctions::cassert(potentialHandle,ISCRITICAL,"atomicSystem::setupPotential: could not open shared library"+libraryFile,__FILE__,__LINE__);
//dlerror();
//verbosity(_Pa,"get the handle to the dualPotential function",2,__FILE__,__LINE__);
//typedef arma::mat (*dualPotential_t)(operatorStruct,paramStruct,mat,cx_mat,mat,mat);
//dualPotential_t dualPotential = (dualPotential_t) dlsym(potentialHandle,"getDualPotential");
//verbosity(_Pa,"check for errors",2,__FILE__,__LINE__);
//const char *dlsym_error = dlerror();
//cout << " the dlsym_error is " << dlsym_error << endl;
//myFunctions::cassert(!dlsym_error,ISCRITICAL,"atomicSystem::setupPotential: Cannot load symbol 'dualPotential'",__FILE__,__LINE__);
////mat *Vdu;
////dualPotential(_Op,_Pa,_dr,_Sf,_X,_G2,Vdu);
////cout << *Vdu << endl;
////cout << dualPotential(_Pa.Z,_dr,_Sf,_X,_G2) << endl;
//verbosity(_Pa,"compute the dual potential using the compiled function ",2,__FILE__,__LINE__);
//mat returnMatrix = dualPotential(_Op,_Pa,_dr,_Sf,_X,_G2);
//verbosity(_Pa,"now write the output matrix "+sourceFile,2,__FILE__,__LINE__);
//this->_Vdual=returnMatrix;
verbosity(_Pa,"atomicSystem::setupPotential: _Vdual has  "+std::to_string(_Vdual.n_rows)+" rows ",2,__FILE__,__LINE__);
verbosity(_Pa,"atomicSystem::setupPotential: _Vdual has  "+std::to_string(_Vdual.n_cols)+" columns ",2,__FILE__,__LINE__);
verbosity(_Pa,"atomicSystem::setupPotential: _Vdual has  "+std::to_string(_Vdual.n_elem)+" elements ",2,__FILE__,__LINE__);
verbosity(_Pa,"atomicSystem::setupPotential: _Vdual has mean "+std::to_string(arma::mean(arma::mean(_Vdual)))+" mean value ",2,__FILE__,__LINE__);
verbosity(_Pa,"atomicSystem::setupPotential: _Vdual has "+std::to_string(arma::max(arma::max(_Vdual)))+" max value ",2,__FILE__,__LINE__);

myFunctions::cassert(((this->_Vdual.is_finite()) && (!this->_Vdual.has_nan())),ISCRITICAL,"atomicSystem::setupPotential: found nan or infinite value in dual potential",__FILE__,__LINE__);
//verbosity(_Pa,"atomicSystem::setupPotential: close dynamic linked library "+sourceFile,2,__FILE__,__LINE__);
//dlclose(potentialHandle);
return 0;

}

int atomicSystem::getAtomicCoordinates()
{

    //! \brief retrive atomic coordinates as specified in .atoms file

//    int returnInt=0;
//
//    verbosity(_Pa,"Now get the position of the atoms in the system and store in matrix",2,__FILE__,__LINE__);
//    mat Xt=this->readAtomicCoordinates();
//    verbosity(_Pa,"We have read "+std::to_string(Xt.n_rows)+" atoms!",2,__FILE__,__LINE__);
//
//    if(_Xt.n_rows==0 || _Xt.n_cols!=3)
//    {
//    returnInt=1;
//    verbosity(_Pa,"Inputfile contains no atomic coordinates - create a sample one",1,__FILE__,__LINE__);
//    this->writeSampleAtomicCoordinateFile();
//    verbosity(_Pa,"Now read that file again as input",1,__FILE__,__LINE__);
//    Xt=this->readAtomicCoordinates();
//    cassert(Xt.n_rows==0,ISCRITICAL,"Sample input file erroneous - no matrix read - contact author or check source code!",__FILE__,__LINE__);
//    verbosity(_Pa,"Now we have a sample atomic system with "+std::to_string(Xt.n_rows)+" atoms!",1,__FILE__,__LINE__);
//    }

    _X=_Pa.X;//!<atomic coordinates stored in columns
    myFunctions::cassert(_X.n_elem>0,ISCRITICAL,"atomicSystem::getAtomicCoordinates - no atoms found in parameter struct",__FILE__,__LINE__);
    return 0;
}

int atomicSystem::setupOptimizers()
{
    //! \brief setup of optimizers for minimizing the nonlinear energy functional
    //! \param Vdual: matrix holds dual potential - lives in real space - enters the electron energy as follows: $E_{ion}=Vdual^{\dagger} n$

    cassert(_Vdual.n_elem>0,ISCRITICAL,"atomicSystem::setupOptimizers potential Vdual has no elements - forgotten to initialize it?",__FILE__,__LINE__);
        verbosity(_Pa,"atomicSystem::setupOptimizers start",2,__FILE__,__LINE__);

    //sdOptimizer sdn(_Op,_Pa,_Vdual,_ltX);
    _sd.reset(new sdOptimizer(_Op,_Pa,_Vdual,_ltX));
            verbosity(_Pa,"atomicSystem::setupOptimizers documentation",2,__FILE__,__LINE__);
    _sd->myLatexClass->commentMyFunction();
    this->_pccg.reset (new pccgOptimizer(_Op,_Pa,_Vdual,_G2,_ltX));
    this->_pccg->myLatexClass->commentMyFunction();
};

int atomicSystem::solveIt()
{
//! \brief minimizing the energy functional by application of steepest descent method first - then conjugate gradient method

cassert(_W->n_rows>0,ISCRITICAL,"atomicSystem::solveIt wavefunction has not been initialized - forgotten to do setupOptimizers!",__FILE__,__LINE__);
cassert(_W->n_cols>0,ISCRITICAL,"atomicSystem::solveIt wavefunction has not been initialized - forgotten to do setupOptimizers!",__FILE__,__LINE__);

_sd->setup(_W);
_sd->myLatexClass->commentMyFunction();
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);
verbosity(_Pa,"atomic: optimize using sd",0,__FILE__,__LINE__);
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);
double finalEnergy=_sd->optimize(_W);
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);
verbosity(_Pa,"atomic: sd optimized system to an energy of "+std::to_string(finalEnergy)+" hartree!",0,__FILE__,__LINE__);
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);
_sd->myLatexClass->commentMyFunction();
verbosity(_Pa,"atomic: restart as orthonormal function",2,__FILE__,__LINE__);
this->orthogonalizeWavefunction();
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);
verbosity(_Pa,"atomic: optimize using pccg",0,__FILE__,__LINE__);
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);
_pccg->myLatexClass->commentMyFunction();
finalEnergy=_pccg->optimize(_W);
_pccg->myLatexClass->commentMyFunction();
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);
//verbosity(_Pa,"atomic: pccg optimized system to an energy of "+std::to_string(finalEnergy)+" hartree!",0,__FILE__,__LINE__);
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);
verbosity(_Pa,"atomic: finished optimization",0,__FILE__,__LINE__);
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);

return 0;
}

atomicSystem::atomicSystem(operatorStruct Op,paramStruct Pa,gnuPlotPlotting *gpLT,latexComment *latX)
{

    //! \brief initialize variables and parser class instance

    this->_caseName=Pa.caseName;
    this->_Pa=Pa;
    this->_Op=Op;
    this->_mat_center_of_cell=arma::ones(_Pa.prodS,3);
    this->_ltX=latX;
    this->_Vdual=arma::zeros(_Pa.prodS,1);

    _gpLT=gpLT;

}

atomicSystem::~atomicSystem()
{
}

int atomicSystem::setupGeometry()
{

//! \brief set up the computational geometry - the grid - see T.Arias' lecture for further explanation

using namespace arma;

string lString="$\\mathbf{S}: 3 \\times 1$ number of sample points along each lattice vector ";
_ltX->newLine(lString);

lString="$ \\Bigg( \\begin{array}{ccc}"+std::to_string(_Pa.S[0])+" \& 0 \& 0 \\\\ ";
lString+=" 0 \& "+ std::to_string(_Pa.S[1]) +" \& 0  \\\\ ";
lString+=" 0 \& 0 \& " + std::to_string(_Pa.S[2]) + "\\end{array} \\Bigg) $ ";
_ltX->newLine(lString);

Col<double> ms(_Pa.prodS);
for(unsigned int i=0;i<ms.n_elem;i++)
{ ms[i]=i; }

Col<double> m1(_Pa.prodS);
Col<double> m2(_Pa.prodS);
Col<double> m3(_Pa.prodS);

for(unsigned int i=0;i<ms.n_elem;i++)
{
m1[i]=std::fmod(ms(i),_Pa.S[0]);
m2[i]=std::floor(std::fmod(ms(i)/_Pa.S[0],_Pa.S[1]));
m3[i]=std::floor(std::fmod(ms(i)/(_Pa.S[0]*_Pa.S[1]),_Pa.S[2]));
}


Mat<double> M(_Pa.prodS,3);
M.fill(0);
M.col(0)=m1;
M.col(1)=m2;
M.col(2)=m3;

Col<double> n1(_Pa.prodS);
Col<double> n2(_Pa.prodS);
Col<double> n3(_Pa.prodS);


for(unsigned int j=0;j<m1.n_rows;j++)
    {
        //value = (expression) ? (if true) : (if false);
        n1[j] = (m1[j]<=(_Pa.S[0]/2)) ? (m1[j]) : (-_Pa.S[0]+m1[j]);
        n2[j] = (m2[j]<=(_Pa.S[1]/2)) ? (m2[j]) : (-_Pa.S[1]+m2[j]);
        n3[j] = (m3[j]<=(_Pa.S[2]/2)) ? (m3[j]) : (-_Pa.S[2]+m3[j]);

        }

Mat<double> N(_Pa.prodS,3);
N.fill(0);
N.col(0)=n1;
N.col(1)=n2;
N.col(2)=n3;

lString=" $\\mathbf{R}: 3 \\times 3$ matrix, each of whose columns is a lattice vector $ R_{k},\\mathbf{R}=\[R_{1},R_{2}, R_{3}\] $ ";
_ltX->newLine(lString);
lString="$ \\Bigg( \\begin{array}{ccc}"+std::to_string(_Pa.R[0,0])+" \& 0 \& 0 \\\\ ";
lString+=" 0 \& "+ std::to_string(_Pa.R[1,1]) +" \& 0  \\\\ ";
lString+=" 0 \& 0 \& " + std::to_string(_Pa.R[2,2]) + "\\end{array} \\Bigg) $ ";
_ltX->newLine(lString);

Mat<double> diagS = diagmat(_Pa.S);
Mat<double> idiagS = inv((mat)diagS);

_ltX->emptyLine();
lString=" $\\mathbf{r}: S_{k} \\times 3$ matrix, matrix, each of whose rows contains the real-space coordinates of a sample point"
" standard order $\\mathbf{r} ≡ \\mathbf{M} Diag(\\vec(S))^{-1}\\mathbf{R}^{T}$  ";
_ltX->newLine(lString);

_r = M*idiagS*_Pa.R.t();

verbosity(_Pa,"compute center of cell",2,__FILE__,__LINE__);

for(int j=0;j<_Pa.prodS;j++)
    {
     _mat_center_of_cell.row(j)=(arma::sum(_Pa.R,1)/2).t();
    }

_dr=arma::sqrt(arma::sum((_r-_mat_center_of_cell)%(_r-_mat_center_of_cell),1));

_ltX->emptyLine();
_ltX->newLine(" The G operator: ");
_ltX->newLine(" components of reciprocal lattice vector: ");
_ltX->newLine(" $\\Pi_{k} S_{k} \\times 3 matrix$ ");
_ltX->newLine(" $ \\mathbf{G} = 2 \\pi \\mathbf{N} \\mathbf{R}^{-1}$ ");

_G=2*arma::datum::pi*N*inv(_Pa.R);

_ltX->emptyLine();
_ltX->newLine(" The G2 operator: ");
_ltX->newLine(" components of square magnitude of the vector in the corresponding row of G: ");
_ltX->newLine(" $\\Pi_{k} S_{k} \\times 1 matrix$ ");
_ltX->newLine(" $ G2 = G^{2}$ ");

_G2=arma::sum(_G%_G,1);

verbosity(_Pa,"now compute the structure factor",2,__FILE__,__LINE__);

verbosity(_Pa,"Compute structure factor using G: "+std::to_string(_G.n_rows)+" x "+std::to_string(_G.n_cols),2,__FILE__,__LINE__);
verbosity(_Pa,"Compute structure factor using X: "+std::to_string(_X.n_rows)+" x "+std::to_string(_X.n_cols),2,__FILE__,__LINE__);

myFunctions::cassert(_X.n_elem>0,ISCRITICAL,"no atomic coordinates read so far - execute readAtomicCoordinates before setupGeometry!",__FILE__,__LINE__);

arma::mat tmpmat(-_G*_X);

verbosity(_Pa,"Now compute exponent of temporary matrix: "+std::to_string(tmpmat.n_rows)+" x "+std::to_string(tmpmat.n_cols),2,__FILE__,__LINE__);

arma::cx_mat exponent(tmpmat.n_rows,tmpmat.n_cols,fill::zeros);
exponent.set_imag(tmpmat);

_ltX->emptyLine();
_ltX->newLine(" The Sf operator: ");
_ltX->newLine(" structure factor associated with the wave vector of the corresponding column of G: ");
_ltX->newLine(" $ \\Pi_{k} S_{k} \\times 1$ column vector (one column matrix) ");
_ltX->newLine(" $ Sf(\\vec(G)) = \\Sigma_{I} \\exp(-i \\vec(G) \\vec(X)_{I}) $ ");

verbosity(_Pa,"Sf is a sum over all atoms I of exp(G*X_{I})!",2,__FILE__,__LINE__);
verbosity(_Pa,"Sf is a column vector Sx1!",2,__FILE__,__LINE__);

_Sf=sum(exp(exponent),1);

verbosity(_Pa,"====================================",0,__FILE__,__LINE__);
verbosity(_Pa,"atomic: setup finished",0,__FILE__,__LINE__);
verbosity(_Pa,"====================================",0,__FILE__,__LINE__);


return 0;
}

int atomicSystem::postDensity(const string texCaption,const string fileNameEnding)
{

//! \brief compute and output electron density
//! \param texCaption: caption of figure in latex document
//! \param fileNameEnding: file ending of output file: e.g. _start or _final - has to be unique in this program

        cx_mat Ui=Uinvers(_Op,_Pa,*_W);
        mat density=computeDensityFromWavefuncs(_Op,_Pa,*_W,Ui);
        //gnuPlotPlotting *gpL;
        string imgName=_gpLT->plotAMatrixSlice(_Pa,_Pa.caseName+"_density_"+fileNameEnding,density,_Pa.S,0);
        _ltX->insertImage(imgName,texCaption);
        //delete gpL;
}

int atomicSystem::postVdual(const string texCaption)
{

//! \brief compute and output Vdual
//! \param texCaption: caption of figure in latex document

        //gnuPlotPlotting *gpL(new gnuPlotPlotting);
        verbosity(_Pa,"atomicSystem::postVdual: plot matrix slice ",2,__FILE__,__LINE__);
        string name=_Pa.caseName+"_Vdual";
        string imgName=_gpLT->plotAMatrixSlice(_Pa,name,this->_Vdual,_Pa.S,0);
        verbosity(_Pa,"atomicSystem::postVdual: now insert image and caption into latex document ",2,__FILE__,__LINE__);
        _ltX->insertImage(imgName,texCaption);
        //delete gpL;
        return 0;
}

int atomicSystem::postPsi()
{

//! \brief compute and output wavefunctions and their squares

        verbosity(_Pa,"atomic: now compute psi",2,__FILE__,__LINE__);
        std::shared_ptr<arma::cx_mat> Psi(new arma::cx_mat(*_W));
        std::shared_ptr<arma::mat> Epsilon(new arma::mat(_Pa.number_of_wavefunctions,1));
        gnuPlotPlotting *gpL;

       if(getPsi(_Op,_Pa,*_W,_Vdual,Psi,Epsilon))
        {
            cout << "succesfully extracted eigenvalues and eigenvectors from the final solution W!" << std::endl;
             cout << "schroedinger: returnEpsilon " << Epsilon->col(0) << endl;

        }
        _ltX->newLine("The electron states have the following energies:");

        arma::mat dat(Psi->n_rows,1);
        string name;
        for(unsigned int st=0; st<Psi->n_cols; st++)
        {
            std::cout << "=== State " << st << ", has energy = " << Epsilon->row(st) << " === " << std::endl;
            _ltX->newLine("=== state " + std::to_string(st) + " has energy " + std::to_string(as_scalar(Epsilon->row(st)))+ "hartree ===");
            verbosity(_Pa,"now we get the square of the wavefunction for output",2,__FILE__,__LINE__);
            dat=real(pow(*_Op.I*Psi->col(st),2));
            verbosity(_Pa,"print three slices of result to as ppm / gnuplot",2,__FILE__,__LINE__);
            #ifdef PLOT_RESULT
            for(int k=0; k<3; k++)
            {
                string plane;
                if(k==0) {plane="yz";}
                else if(k==1) {plane="xz";}
                else if(k==2) {plane="xy";}
                verbosity(_Pa,"now take a slice of data at plane "+plane,2,__FILE__,__LINE__);
                arma::mat sl;
                verbosity(_Pa,"take a slice "+plane,2,__FILE__,__LINE__);
                sl=myFunctions::slice(_Pa,dat,_Pa.S,(int) _Pa.S(k)/2.,k);
                #ifdef PLOT_PPM
                verbosity(_Pa,"get name for ppm file "+plane,2,__FILE__,__LINE__);
                name="psi"+std::to_string(st)+"d_m_"+plane; //std::to_string(k)
                verbosity(_Pa,"write ppt file "+name,2,__FILE__,__LINE__);
                ppm(name,sl*0.3,sl,sl,&latX);
                #endif
                #ifdef PLOT_GNUPLOT
                if(st==0) _ltX->newLine("The physical eingenstates of the hamiltonian Psi are shown in the next figure:");
                if(k==0)
                {
                    name="schroedinger_psi_"+std::to_string(st);
                    string imageName=_gpLT->plotAMatrixSlice(_Pa,name,dat,_Pa.S,k);
                    string caption="result from solution of schroedinger equation: slice through wave function of electron " + std::to_string(st) + " plane "+plane;
                    _ltX->insertImage(imageName,caption);
                }
                #endif
            }
            #endif
        }
        return 0;
}

int atomicSystem::setupWavefunction()
{
//! \brief generate inital wavefunction as random values

    _W.reset (new arma::cx_mat(arma::randn<arma::cx_mat>(_Pa.prodS,_Pa.number_of_wavefunctions)));

}