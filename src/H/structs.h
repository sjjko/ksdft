#ifndef STRUCTS_H
#define STRUCTS_H

#include "O.h"
#include "Linv.h"
#include "L.h"
#include "cJ.h"
#include "cI.h"
#include "cJdag.h"
#include "cIdag.h"
#include "main.h"

using namespace arma;

// some forward declarations
class cI;
class cIdag;
class cJ;
class cJdag;
class Lop;
class Linv;
class OOp;

struct paramStruct
{
    vec S; //!< the dimensional vector
    mat R;
    mat X; //!< atomic coordinates
    double Z;  //!< charge number of atom
    double f;  //!< number of electrons per atom
    int number_of_wavefunctions;  //!< number of wavefunction per electron
    double alpha;  //!< step parameter for steepest descent method
    int sdNit; //!< number of iterations steepest descent
    int pccgNit; //!< number of iterations conjugate gradient
    int globalVL; //!< global verbosity level
    double prodS;
    string caseName;
};

struct operatorStruct
{
    int testint2;
    std::shared_ptr<int> testint;
    std::shared_ptr<OOp> O;
    std::shared_ptr<Lop> L;
    std::shared_ptr<Linv> Li;
    std::shared_ptr<cJ> J;
    std::shared_ptr<cJdag> Jd;
    std::shared_ptr<cI> I;
    std::shared_ptr<cIdag> Id;
};



#endif
