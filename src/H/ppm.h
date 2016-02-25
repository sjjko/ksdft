#ifndef PPM_H
#define PPM_H

/*
 * =====================================================================================
 *
 *       Filename:  ppm.h
 *
 *    Description:  ppm.h
 *
 *        Version:  1.0
 *        Created:  2015-03-09 00:33:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  T.Arias - converted to c++ by S.KonzettStoffl
 *   Organization:
 *
 * =====================================================================================
 */


#include "main.h"
#include "smooth.h"

using namespace std;
using namespace myFunctions;

/*  Usage: ppm(fname,red,green,blue)
%#
%# Output- color ppm image in file "fname" (view with "xli fname").
%# Input- red, green, blue: 2d data of red, green, blue intensities
*/
inline int ppm(string fname,arma::mat red,arma::mat green,arma::mat blue,latexComment *latX)
{
    //! \brief output routine for ppm images
    //!
    //! \param red, green, blue: 2d data of red, green, blue intensities
    //! \param
    //! \return color ppm image in file "fname" (view with "xli fname").


ofstream myfile;
myfile.open (fname+".ppm", std::ofstream::out);
//%# Enlarge image
verbosity("write ppt file ",2,__FILE__,__LINE__);

arma::mat redSmooth,greenSmooth,blueSmooth;

//for (int en=1;en<=4;en++)
//{
    verbosity("smooth the datasets using the smooth routine ",2,__FILE__,__LINE__);
    redSmooth=smooth(red);
    greenSmooth=smooth(green);
    blueSmooth=smooth(blue);
//}
    verbosity("ppm: finished smoothing ",2,__FILE__,__LINE__);

cassert(!red.has_inf(),ISCRITICAL,"red has inf",__FILE__,__LINE__);
cassert(!red.has_nan(),ISCRITICAL,"red has nan",__FILE__,__LINE__);
cassert(red.is_finite(),ISCRITICAL,"red has infinite values",__FILE__,__LINE__);

double pixmx=255;
double height=size(redSmooth,0);
double width=size(redSmooth,1);
    verbosity("ppm: output has size of smoothed datamatrix size of height "+std::to_string(height),2,__FILE__,__LINE__);

//arma::vec<arma::mat<double>> tmpv(red,green,blue);
double mx=0;
verbosity("ppm: get maximal values ",2,__FILE__,__LINE__);

mx=max(mx,redSmooth.max());
mx=max(mx,blueSmooth.max());
mx=max(mx,greenSmooth.max());
double mn=1000000000;
mn=min(mn,redSmooth.min());
mn=min(mn,blueSmooth.min());
mn=min(mn,greenSmooth.min());
//fid=fopen(fname,’w’);
myfile << "P3" << endl;
//fprintf(fid,’P3\n’);
myfile << width << " " << height << endl;
//fprintf(fid,’%d %d\n’,width,height);
myfile << pixmx << endl;
//fprintf(fid,’%d\n’,pixmx);

verbosity("ppm: reshape colour arrays ",2,__FILE__,__LINE__);
arma::Col<double> colRed(reshape(redSmooth,redSmooth.n_rows*redSmooth.n_cols,1));
arma::Col<double> colGreen(reshape(greenSmooth,greenSmooth.n_rows*greenSmooth.n_cols,1));
arma::Col<double> colBlue(reshape(blueSmooth,blueSmooth.n_rows*blueSmooth.n_cols,1));

verbosity("ppm: now fill into data columns ",2,__FILE__,__LINE__);

arma::mat ppmDat(colBlue.n_elem,3);
ppmDat.col(0)=colRed;
ppmDat.col(1)=colGreen;
ppmDat.col(2)=colBlue;
arma::mat ones=ppmDat.ones();

verbosity("ppm: now rescale the data using maximum mx and minimum of data mn",2,__FILE__,__LINE__);
ppmDat=(ppmDat-ones*mn)/(ones*mx-ones*mn)*pixmx; //!< now rescale the data according to min/max values
verbosity("ppm: round the values",2,__FILE__,__LINE__);
//dat=([arma::reshape(red’,1,width*height); ...
//		arma::reshape(green’,1,width*height); ...
//		arma::reshape(blue’,1,width*height)] ...
//		-mn)/(mx-mn)*pixmx;
verbosity("ppm: output the scaled data",2,__FILE__,__LINE__);

cassert(!ppmDat.has_inf(),ISCRITICAL,"ppmDat has inf",__FILE__,__LINE__);
cassert(!ppmDat.has_nan(),ISCRITICAL,"ppmDat has nan",__FILE__,__LINE__);
cassert(ppmDat.is_finite(),ISCRITICAL,"ppmDat has infinite values",__FILE__,__LINE__);

//myfile.setf(ios::fixed,ios::right);
myfile << ppmDat << std::endl;
//fprintf(fid,’%d ’,dat);
verbosity("ppm: close the data file",2,__FILE__,__LINE__);

myfile.close();

#ifdef INCLUDE_PICTURE_IN_LATEX
latX->latexComment::convertPPMToPS(fname);
latX->insertImage(fname+".ps"); //!< now we insert the image into our latex document
#endif

return 1;
}


#endif // PPM_H
