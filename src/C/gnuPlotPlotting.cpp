#include "gnuPlotPlotting.h"

using namespace myFunctions;

gnuPlotPlotting::gnuPlotPlotting()
{
    //ctor
}

gnuPlotPlotting::~gnuPlotPlotting()
{
    //dtor
}

int gnuPlotPlotting::plotMatrix3Slices(string title,arma::mat matrixToPlot,Col<double> S)
{

  //! \brief does not plot - takes only slices and gives output in .2dMat files which are automatially postprocesses by
  //! script plotScript.sh
    //!
    //! \param title: title string for matrix
    //! \param matrixToPlot real arma matrix for slicing and output -> Sx1 data!
    //! \return nothing (int=0)

arma::mat TmpMat;
for(int k=0;k<3;k++)
{
    string plane;
    if(k==0) {plane="yz";}
    else if(k==1) {plane="xz";}
    else if(k==2) {plane="xy";}
    TmpMat=myFunctions::slice(matrixToPlot,S,S(k)/2.,k);
    string fileName;
    fileName = title +"_"+ plane +".2dMat";
    ofstream of(fileName);
    of << TmpMat << endl;
    of.close();
}

return 0;
}


std::string gnuPlotPlotting::plotAMatrixSlice(string title,arma::mat matrixToPlot,Col<double> S,int sliceIndex)
{

  //! \brief does not plot - takes only slices and gives output in .2dMat files which are automatially postprocesses by
  //! script plotScript.sh
    //!
    //! \param title: title string for matrix
    //! \param matrixToPlot real arma matrix for slicing and output -> Sx1 data!
    //! \return nothing (int=0)

    arma::mat TmpMat;

    if(sliceIndex==0) {_nameOfActualPlane="yz";}
    else if(sliceIndex==1) {_nameOfActualPlane="xz";}
    else if(sliceIndex==2) {_nameOfActualPlane="xy";}
    TmpMat=myFunctions::slice(matrixToPlot,S,S(sliceIndex)/2.,sliceIndex);
    string fileName;
    fileName = title +"_"+ _nameOfActualPlane +".2dMat";

    ofstream of(fileName);
    of << TmpMat << endl;
    of.close();

    return title +"_"+ _nameOfActualPlane;
}

