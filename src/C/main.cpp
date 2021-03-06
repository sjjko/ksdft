/** density functional theory code
 * devised after T. Arias lecture on dft for matlab
 * by Stefan Konzett Stoffl
 * 2015-2016

 for further info see http://sjjko.github.io/
 */

#include "main.h"
#ifdef WITH_TEX
 #ifdef USE_EXTERNAL_LIB_TEXCALLER
    #include <texcaller.h>
#endif // USE_EXTERNAL_LIB_TEXCALLER
    //include to generate the image folder structure
    #include <sys/types.h>
    #include <sys/stat.h>
#endif
#include "latexComment.h"
#include "structs.h"
#include "O.h"
#include "L.h"
#include "Linv.h"
#include "cI.h"
#include "cJ.h"
#include "cJdag.h"
#include "cIdag.h"
#include "structs.h"
#include "optimizerBase.h"
#include "sdOptimizer.h"
#include "pccgOptimizer.h"
#include "smooth.h"
#include "ppm.h"
#include <typeinfo>
#include "inputParser.h"
#include "JItest.h"
#include "slice.h"
#include "testingClass.h"
#include "gnuPlotPlotting.h"
#include "atomicSystem.h"

using namespace std;
using namespace arma;


int main(int argc, char** argv)
{

paramStruct Pa;
timerClass tC=timerClass(Pa);
  boost::timer::auto_cpu_timer t;

#include "inputArguments.h"

#ifdef PLOT_GNUPLOT
    gnuPlotPlotting gpL;
#endif // PLOT_GNUPLOT

#ifdef WITH_TEX
    latexComment latX;
    const std::string latexFileName="ksdft++.pdf";
    latX.writeHeader();
    latX.section(" ksdft++ documentation: ");
    latX.newLine(" Program for computing the electronic energy configuration of atoms and molecules. ");
    latX.newLine(" Devised after lecture by T.Arias. This run computes the case "+Pa.caseName);

    latX.subSection("The inital setup stage");
    latX.newLine(" The following steps are performed in the setup stage: ");
#endif


arma::mat Xt; //!< the atomic coordinate matrix
verbosity(Pa,"Read Input Parameters!",2,__FILE__,__LINE__);
#ifdef WITH_BOOST
    verbosity(Pa,"main: parse input file "+paramFile,2,__FILE__,__LINE__);
    inputParser Iparser(pathToCase,paramFile,&Pa);
    if(Iparser.parseInput()){cout << "successfully read input parameters" << endl;};

    verbosity(Pa,"Now get the position of the atoms in the system and store in matrix",2,__FILE__,__LINE__);
    Xt=Iparser.readAtomicCoordinates();
    verbosity(Pa,"We have read "+std::to_string(Xt.n_rows)+" atoms!",2,__FILE__,__LINE__);

    if(Xt.n_rows==0 || Xt.n_cols!=3)
    {
    verbosity(Pa,"Inputfile contains no atomic coordinates - create a sample one",1,__FILE__,__LINE__);
    Iparser.writeSampleAtomicCoordinateFile();
    verbosity(Pa,"Now read that file again as input",1,__FILE__,__LINE__);
    Xt=Iparser.readAtomicCoordinates();
    cassert(Xt.n_rows==0,ISCRITICAL,"Sample input file erroneous - no matrix read - contact author or check source code!",__FILE__,__LINE__);
    verbosity(Pa,"Now we have a sample atomic system with "+std::to_string(Xt.n_rows)+" atoms!",1,__FILE__,__LINE__);
    }
    arma::mat X=Xt.t();
    Pa.X=X;

#else

    myFunctions::cassert(false,ISCRITICAL,"sorry - not yet implemented - runs only with boost lib!",__FILE__,__LINE__);

//TODO: dysfunctional !!! not yet updated
    paramStruct Pa;
    Pa={.Z=1,.f=2,.number_of_wavefunctions=2,.alpha=0.5e-3,.Nit=2,.globalVL=333};
#endif



verbosity(Pa,"For reasons of practicality we read X rowwise, but compute now columnwise - take transpose of Xt to get X! ",1,__FILE__,__LINE__);

#ifdef WITH_TEX
    latX.emptyLine();
    string lString=" $\\mathbf{X}:$ atomic locations, $N \\times 3$ matrix, each of whose rows contains the real-space coordinates of an atom (in bohrs)";
    latX.newLine(lString);
#endif // WITH_TEX

verbosity(Pa,"start minidft++",2,__FILE__,__LINE__);

verbosity(Pa,"do the setup first!",2,__FILE__,__LINE__);
#include "setup.h"
verbosity(Pa,"setup has been done!",2,__FILE__,__LINE__);

#ifdef WITH_TEX
latX.subSection("Initialization of the operator classes");
#endif // WITH_TEX

std::shared_ptr<OOp> Oclass(new OOp(Pa.R,&latX));
#ifdef WITH_TEX
    latX.startItemize();
    Oclass.get()->myLatexClass->commentMyFunctionAsItem();
#endif
std::shared_ptr<Lop> Lclass(new Lop(G2,G2comp,Pa.R,Pa.number_of_wavefunctions,&latX));
Lclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<Linv> Liclass(new Linv(G2,Pa.R,Pa.number_of_wavefunctions,&latX));
Liclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<cJ> cJclass(new cJ(Pa.S,Pa.number_of_wavefunctions,&latX));
cJclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<cI> cIclass(new cI(Pa,Pa.S,Pa.number_of_wavefunctions,&latX));
cIclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<cIdag> cIdagclass(new cIdag(Pa,Pa.S,Pa.number_of_wavefunctions,&latX));
cIdagclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<cJdag> cJdagclass(new cJdag(Pa.S,Pa.number_of_wavefunctions,&latX));
cJdagclass.get()->myLatexClass->commentMyFunctionAsItem();
#ifdef WITH_TEX
latX.endItemize();
#endif

verbosity(Pa,"now we initiate OOp operator class!",2,__FILE__,__LINE__);
operatorStruct Op;
verbosity(Pa,"initiate OOp!",2,__FILE__,__LINE__);
Op.O=Oclass;
verbosity(Pa,"initiate LOp!",2,__FILE__,__LINE__);
verbosity(Pa,"initiate LinvOp!",2,__FILE__,__LINE__);
Op.L=Lclass;
Op.Li=Liclass;
verbosity(Pa,"initiate cIOp!",2,__FILE__,__LINE__);
Op.I=cIclass;
verbosity(Pa,"initiate IdagOp!",2,__FILE__,__LINE__);
Op.Id=cIdagclass;
verbosity(Pa,"initiate cJOp!",2,__FILE__,__LINE__);
Op.J=cJclass;
verbosity(Pa,"initiate cJdagOp!",2,__FILE__,__LINE__);
Op.Jd=cJdagclass;
verbosity(Pa,"operator struct operated!",2,__FILE__,__LINE__);


#ifdef PERFORM_OPERATOR_TESTS
verbosity(Pa,"We test the implementation of the operators!",0,__FILE__,__LINE__);
testingClass* ts=new testingClass(&latX,Pa);
ts->testForIJ(Op,Pa,Pa.S);
ts->testForIJMultiColumn(Op,Pa,Pa.S);
ts->testPoisson(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell);
ts->testEwald(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell,Sf,X);
ts->testHermitian(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell,Sf,X);
ts->testMultiColumn(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell,Sf,X);
ts->testDiagOuter();
ts->testSchroedinger(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell,Sf,X,G2comp);
verbosity(Pa,"Finished testing - dump the class instance!",0,__FILE__,__LINE__);
delete ts;
#endif

verbosity(Pa,"Main computing part:",1,__FILE__,__LINE__);

atomicSystem atm(Op,Pa,&gpL,&latX);
verbosity(Pa,"main: retrieve atomic coordinates",1,__FILE__,__LINE__);
atm.getAtomicCoordinates();
atm.setupGeometry();
atm.setupPotential();
atm.postVdual("dual ion potential");
verbosity(Pa,"main: setup wavefunction",1,__FILE__,__LINE__);
atm.setupWavefunction();
verbosity(Pa,"main: setup optimizers",1,__FILE__,__LINE__);
atm.setupOptimizers();
verbosity(Pa,"main: solve and retrieve the final Solution per string",1,__FILE__,__LINE__);

string finalSolution = atm.solveIt();
verbosity(Pa,"main: do postprocessing",1,__FILE__,__LINE__);
finalSolution += atm.doPostProcessing("final density after mininimizing the energy functional","final");
verbosity(Pa,"main: output final energy result in text file",1,__FILE__,__LINE__);
string finalResultFilename="finalResult_"+Pa.caseName+".dat";
std::ofstream finalResultFile(finalResultFilename);
finalResultFile << finalSolution << endl;
finalResultFile.close();

verbosity(Pa,"Latex part:",2,__FILE__,__LINE__);

#ifdef WITH_TEX
    std::ofstream latexStringFile("latex.tex");
    latX.closeString();
    latexStringFile << latX.getString() << endl;
    latexStringFile.close();
#endif // WITH_TEX

std::cout << "ksdft++ is finished!" << std::endl;

return(0);
}
