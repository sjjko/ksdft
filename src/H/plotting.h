#ifndef PLOTTING_H
#define PLOTTING_H

#include "armadillo"
#include "plotImplementationBase.h"


class plotting
{
    public:
        plotting()
        {
        };
        virtual ~plotting(){};
        int plotMatrix(std::string title,arma::mat matrixToPlot);
        plotImplementationBase plot();
       // typedef typename plotImplementationBase::plot() plot();

};

#endif // PLOTTING_H
