#include "openCVPlotting.h"

using namespace cv;
using namespace myFunctions;

openCVPlotting::openCVPlotting()
{
}

openCVPlotting::~openCVPlotting()
{
    //dtor
}

int openCVPlotting::plotMatrix(string title,arma::mat matrixToPlot)
{

    verbosity("openCVPlotting: matrix row order in cv and column order in armadillo - have to transpose first",1,__FILE__,__LINE__);
    arma::mat matrixToPlott=matrixToPlot.t();
    cv::Mat opencv_mat( matrixToPlott.n_rows, matrixToPlott.n_cols, CV_64FC1, matrixToPlott.memptr() );
    cv::Mat image;

    cassert(!matrixToPlott.has_inf(),ISCRITICAL,"matrixToPlott has inf",__FILE__,__LINE__);
    cassert(!matrixToPlott.has_nan(),ISCRITICAL,"matrixToPlott has nan",__FILE__,__LINE__);

    cv::applyColorMap( opencv_mat, image, cv::COLORMAP_HOT );
    cv::imshow("colormap", image);
    cv::waitKey(0);

}
