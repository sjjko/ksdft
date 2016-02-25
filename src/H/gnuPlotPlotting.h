#ifndef GNUPLOTPLOTTING_H
#define GNUPLOTPLOTTING_H

#include "main.h"
#include "customAssert.h"
#include "slice.h"

using namespace std;

class gnuPlotPlotting
{
    private:
        string _nameOfActualPlane;

    public:

        gnuPlotPlotting();
        virtual ~gnuPlotPlotting();
        int plotMatrix3Slices(string title,arma::mat matrixToPlot,Col<double> S);
        string plotAMatrixSlice(string title,arma::mat matrixToPlot,Col<double> S,int sliceIndex);
        inline string getPlaneName(){return _nameOfActualPlane;};
};

#endif // GNUPLOTPLOTTING_H
