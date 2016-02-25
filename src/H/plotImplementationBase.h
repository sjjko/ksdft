#ifndef PLOTIMPLEMENTATIONBASE_H
#define PLOTIMPLEMENTATIONBASE_H


#include "armadillo"
#include "customAssert.h"


class plotImplementationBase
{
    public:
        plotImplementationBase(){};
        virtual ~plotImplementationBase(){};
        virtual int plotMatrix(std::string title,arma::mat matrixToPlot){};

};

#endif // PLOTIMPLEMENTATIONBASE_H
