#include "plotting.h"

//template<class plotImplementationBase> plotting<plotImplementationBase>::plotting()
//{
//    //ctor
//}
//
//template<class plotImplementationBase> plotting<plotImplementationBase>::~plotting()
//{
//    //dtor
//}

int plotting<plotImplementationBase>::plotMatrix(string title,arma::mat matrixToPlot)
{

plotImplementationBase pL = new plotImplementationBase;//(title,matrixToPlot);
pL.plot();

return 0;

}
