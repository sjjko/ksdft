/** \brief density functional theory code
 * devised after T. Arias lecture on dft for matlab
 * by Stefan Konzett Stoffl
 * 02-2016
 */
#include "main.h"
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
#include "structs.h"
#include "inputParser.h"
#include "latexComment.h"
#include "JItest.h"
#include "slice.h"
#include "testingClass.h"
#include "gnuPlotPlotting.h"

using namespace std;
using namespace arma;

int GLOBAL_verbosityLevel=0;

int main(int argc, char** argv)
{

verbosity("start with "+std::to_string(argc)+"input arguments.",2,__FILE__,__LINE__);
myFunctions::cassert(argc==2,ISCRITICAL,"give case directory as first argument to ksdft++!",__FILE__,__LINE__);
string pathToCase=argv[1];
verbosity("check if we have a true case directory at "+pathToCase,2,__FILE__,__LINE__);
string paramFile=pathToCase+"/inputFile.param";
ifstream f(paramFile.c_str());
myFunctions::cassert(f.good(),ISCRITICAL,"could not find a inputFile.param file in above directory - Exit!",__FILE__,__LINE__);
f.close();
verbosity("our input file for this run is "+paramFile,2,__FILE__,__LINE__);

#ifdef PLOT_GNUPLOT
    gnuPlotPlotting gpL;
#endif // PLOT_GNUPLOT

#ifdef WITH_TEX
    latexComment latX;
    const std::string latexFileName="ksdft++.pdf";
    latX.writeHeader();
    latX.section(" ksdft++ documentation: ");
    latX.subSection("The inital setup stage");
    latX.newLine(" The following steps are performed in the setup stage: ");
#endif

//registryClass *Reg=new registryClass();
//errorClass *ErrO=new errorClass(Reg);
//checkOperatorSize<arma::vec> *checkVec = new checkOperatorSize<arma::vec>();

arma::mat Xt; //!< the atomic coordinate matrix
verbosity("Read Input Parameters!",2,__FILE__,__LINE__);
#ifdef WITH_BOOST
    paramStruct Pa;
    inputParser Iparser(paramFile,&Pa);
    if(Iparser.parseInput()){cout << "successfully read input parameters" << endl;};

    verbosity("Now get the position of the atoms in the system and store in matrix",2,__FILE__,__LINE__);
    Xt=Iparser.readAtomicCoordinates();
    verbosity("We have read "+std::to_string(Xt.n_rows)+" atoms!",2,__FILE__,__LINE__);

    if(Xt.n_rows==0 || Xt.n_cols!=3)
    {
    verbosity("Inputfile contains no atomic coordinates - create a sample one",1,__FILE__,__LINE__);
    Iparser.writeSampleAtomicCoordinateFile();
    verbosity("Now read that file again as input",1,__FILE__,__LINE__);
    Xt=Iparser.readAtomicCoordinates();
    cassert(Xt.n_rows==0,ISCRITICAL,"Sample input file erroneous - no matrix read - contact author or check source code!",__FILE__,__LINE__);
    verbosity("Now we have a sample atomic system with "+std::to_string(Xt.n_rows)+" atoms!",1,__FILE__,__LINE__);
    }


#else
//TODO: dysfunctional !!! not yet updated
    paramStruct Pa;
    Pa={.Z=1,.f=2,.number_of_wavefunctions=2,.alpha=0.5e-3,.Nit=2,.globalVL=333};
#endif

verbosity("For reasons of practicality we read X rowwise, but compute now columnwise - take transpose of Xt to get X! ",1,__FILE__,__LINE__);

#ifdef WITH_TEX
    arma::mat X=Xt.t();
    latX.emptyLine();
    string lString=" $\\mathbf{X}:$ atomic locations, $N \\times 3$ matrix, each of whose rows contains the real-space coordinates of an atom (in bohrs)";
    latX.newLine(lString);
#endif // WITH_TEX

verbosity("start minidft++",2,__FILE__,__LINE__);

verbosity("do the setup first!",2,__FILE__,__LINE__);
#include "setup.h"
verbosity("setup has been done!",2,__FILE__,__LINE__);

#ifdef WITH_TEX
latX.subSection("Initialization of the operator classes");
#endif // WITH_TEX

std::shared_ptr<OOp> Oclass(new OOp(Pa.R,&latX));
#ifdef WITH_TEX
    latX.startItemize();
    Oclass.get()->myLatexClass->commentMyFunctionAsItem();
#endif
std::shared_ptr<Lop> Lclass(new Lop(G2,Pa.R,Pa.number_of_wavefunctions,&latX));
Lclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<Linv> Liclass(new Linv(G2,Pa.R,Pa.number_of_wavefunctions,&latX));
Liclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<cJ> cJclass(new cJ(Pa.S,Pa.number_of_wavefunctions,&latX));
cJclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<cI> cIclass(new cI(Pa.S,Pa.number_of_wavefunctions,&latX));
cIclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<cIdag> cIdagclass(new cIdag(Pa.S,Pa.number_of_wavefunctions,&latX));
cIdagclass.get()->myLatexClass->commentMyFunctionAsItem();
std::shared_ptr<cJdag> cJdagclass(new cJdag(Pa.S,Pa.number_of_wavefunctions,&latX));
cJdagclass.get()->myLatexClass->commentMyFunctionAsItem();
#ifdef WITH_TEX
latX.endItemize();
#endif

verbosity("now we initiate OOp operator class!",2,__FILE__,__LINE__);
operatorStruct Op;
verbosity("initiate OOp!",2,__FILE__,__LINE__);
Op.O=Oclass;
verbosity("initiate LOp!",2,__FILE__,__LINE__);
Op.L=Lclass;
verbosity("initiate LinvOp!",2,__FILE__,__LINE__);
Op.Li=Liclass;
verbosity("initiate cIOp!",2,__FILE__,__LINE__);
Op.I=cIclass;//(S,Pa.number_of_wavefunctions);
verbosity("initiate IdagOp!",2,__FILE__,__LINE__);
Op.Id=cIdagclass; //(S,Pa.number_of_wavefunctions);
verbosity("initiate cJOp!",2,__FILE__,__LINE__);
Op.J=cJclass;//(S,Pa.number_of_wavefunctions);
verbosity("initiate cJdagOp!",2,__FILE__,__LINE__);
Op.Jd=cJdagclass; //(S,Pa.number_of_wavefunctions);
verbosity("operator struct operated!",2,__FILE__,__LINE__);

#ifdef PERFORM_OPERATOR_TESTS
verbosity("We test the implementation of the operators!",0,__FILE__,__LINE__);
testingClass* ts=new testingClass(&latX);
ts->testForIJ(Op,Pa,Pa.S);
ts->testForIJMultiColumn(Op,Pa,Pa.S);
ts->testPoisson(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell);
ts->testEwald(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell,Sf,X);
ts->testHermitian(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell,Sf,X);
ts->testMultiColumn(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell,Sf,X);
ts->testDiagOuter();
ts->testSchroedinger(Op,Pa,Pa.S,Pa.R,r,mat_center_of_cell,Sf,X,G2);
verbosity("Finished testing - dump the class instance!",0,__FILE__,__LINE__);
delete ts;
#endif

#ifdef FDTEST
    #include "fdtest.h"
#endif

//#ifdef SCHROEDINGER
//    #include "schroedinger.h"
//#endif

#ifdef HATOMS
    #include "Hatoms.h"
#endif

#ifdef WITH_TEX
    std::ofstream latexStringFile("latex.tex");
    latX.closeString();
    latexStringFile << latX.getString() << endl;
    latexStringFile.close();
#endif // WITH_TEX

//#ifdef WITH_TEX
//    latX.closeString();
//    latX.writeTheLatexDocument(latexFileName);
//#endif

std::cout << "minidft++ is finished!" << std::endl;

return(0);
}
