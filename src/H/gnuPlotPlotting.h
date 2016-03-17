#ifndef GNUPLOTPLOTTING_H
#define GNUPLOTPLOTTING_H


/*! formatted output of matrix data for postprocessing with plotScript.sh (using gnuplot)
*/

#include "main.h"
#include "customAssert.h"
#include "slice.h"
#include "structs.h"
#include "smooth.h"

using namespace std;

class gnuPlotPlotting
{
    private:
        string _nameOfActualPlane;

    public:

        gnuPlotPlotting();
        virtual ~gnuPlotPlotting();
        int plotMatrix3Slices(paramStruct Pa,string title,arma::mat matrixToPlot,Col<double> S);
        string plotAMatrixSlice(paramStruct Pa,string title,arma::mat matrixToPlot,Col<double> S,int sliceIndex);
        inline string getPlaneName(){return this->_nameOfActualPlane;};
        inline string setPlaneName(string in){_nameOfActualPlane=in;};
};

#endif // GNUPLOTPLOTTING_H
