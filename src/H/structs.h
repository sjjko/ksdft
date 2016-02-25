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


struct paramStruct
{
    vec S; //!< the dimensional vector
    mat R;
    double Z;  //!< charge number of atom
    double f;  //!< number of electrons per atom
    int number_of_wavefunctions;  //!< number of wavefunction per electron
    double alpha;  //!< step parameter for steepest descent method
    int sdNit; //!< number of iterations steepest descent
    int pccgNit; //!< number of iterations conjugate gradient
    int globalVL; //!< global verbosity level
    double prodS;
};

struct operatorStruct
{
    int testint2;
    //OOp OO(arma::mat Ri);
    std::shared_ptr<int> testint;
    std::shared_ptr<OOp> O;//(new OOp(arma::mat Ri));
    std::shared_ptr<Lop> L;
    std::shared_ptr<Linv> Li;
    std::shared_ptr<cJ> J;
    std::shared_ptr<cJdag> Jd;
    std::shared_ptr<cI> I;
    std::shared_ptr<cIdag> Id;
    /*std::shared_ptr<Lop>  L=new Lop(arma::mat G2inp,arma::mat Rinp,int Numwavfunc);
    std::shared_ptr<Linv>  Li=new Linv(arma::mat G2inp,arma::mat Rinp,int Numwavfunc);
    std::shared_ptr<cJ>  J=new cJ(arma::Col<double> Sinp,int Numwavfunc);
    std::shared_ptr<cJdag>  Jd=new cJdag(arma::Col<double> Sinp,int Numwavfunc);
    std::shared_ptr<cI>  I=new cI(arma::Col<double> Sinp,int Numwavfunc);
    std::shared_ptr<cIdag> Id=new cIdag(arma::Col<double> Sinp,int Numwavfunc);*/
};



#endif
