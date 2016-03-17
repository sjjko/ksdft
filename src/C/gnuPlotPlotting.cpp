#include "gnuPlotPlotting.h"

using namespace myFunctions;

gnuPlotPlotting::gnuPlotPlotting()
{
_nameOfActualPlane="";
}

gnuPlotPlotting::~gnuPlotPlotting()
{
    //dtor
}

int gnuPlotPlotting::plotMatrix3Slices(paramStruct Pa, string title,arma::mat matrixToPlot,Col<double> S)
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
    TmpMat=myFunctions::slice(Pa,matrixToPlot,S,S(k)/2.,k);
    string fileName;
    fileName = title +"_"+ plane +".2dMat";
    ofstream of(fileName);
    of << TmpMat << endl;
    of.close();
}

return 0;
}


std::string gnuPlotPlotting::plotAMatrixSlice(paramStruct Pa,string title,arma::mat matrixToPlot,Col<double> S,int sliceIndex)
{

  //! \brief does not plot - takes only slices and gives output in .2dMat files which are automatially postprocesses by
  //! script plotScript.sh
    //!
    //! \param title: title string for matrix
    //! \param matrixToPlot real arma matrix for slicing and output -> Sx1 data!
    //! \return nothing (int=0)

    arma::mat TmpMat;
    myFunctions::cassert(sliceIndex<3 && sliceIndex>-1,ISCRITICAL,"gnuPlotPlotting::plotAMatrixSlice: sliceIndex out of range",__FILE__,__LINE__);
    verbosity(Pa,"gnuPlotPlotting::plotAMatrixSlice set sliceIndex plane"+std::to_string(sliceIndex),2,__FILE__,__LINE__);
    if(sliceIndex==0)
    {
    _nameOfActualPlane="yz";
    }
    else if(sliceIndex==1) {_nameOfActualPlane="xz";}
    else if(sliceIndex==2) {_nameOfActualPlane="xy";}
    verbosity(Pa,"gnuPlotPlotting::plotAMatrixSlice set sliceIndex plane finished",2,__FILE__,__LINE__);
    verbosity(Pa,"gnuPlotPlotting::plotAMatrixSlice now slice the matrix",2,__FILE__,__LINE__);
    TmpMat=myFunctions::slice(Pa,matrixToPlot,S,S(sliceIndex)/2.,sliceIndex);
    string fileName;
    fileName = title +"_"+ _nameOfActualPlane +".2dMat";
    verbosity(Pa,"gnuPlotPlotting::now output into file stream",2,__FILE__,__LINE__);

    verbosity(Pa,"gnuPlotPlotting:: plotAMatrixSlice now smooth the output data",2,__FILE__,__LINE__);


    //mat smoothedMat(smooth(Pa,TmpMat));
    mat TmpIntermediateMat(TmpMat);
    for(int sI=0;sI<Pa.smoothingIterations;sI++)
    {
    mat smoothedTmpMat(smooth(Pa,TmpIntermediateMat));
    TmpIntermediateMat=smoothedTmpMat;
    }

    ofstream of(fileName);
    of << TmpIntermediateMat << endl;
    of.close();

    string returnString=title +"_"+ _nameOfActualPlane;
    return returnString;
}

