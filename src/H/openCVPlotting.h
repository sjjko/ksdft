#ifndef OPENCVPLOTTING_H
#define OPENCVPLOTTING_H

#include "armadillo"
#include "customAssert.h"
#include "plotImplementationBase.h"
#include <opencv2/contrib/contrib.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

class openCVPlotting: public plotImplementationBase
{
    public:
        openCVPlotting();
        virtual ~openCVPlotting();
        int plotMatrix(std::string title,arma::mat matrixToPlot);

};

#endif // OPENCVPLOTTING_H
