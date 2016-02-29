#ifndef SMOOTH_H
#define SMOOTH_H

#include "main.h"
#include "structs.h"

arma::mat smooth(arma::mat in);
arma::mat colsmooth(arma::mat dat);

//%# Double row density by smoothing columns
using namespace myFunctions;

inline arma::mat colsmooth(paramStruct Pa, arma::mat dat)
{
    int nc=dat.n_rows;
    verbosity(Pa,"colsmooth: allocate matrix, which has two times as many rows as the input matrix "+std::to_string(dat.size()),2,__FILE__,__LINE__);
    arma::mat out=arma::zeros<arma::mat>(2*nc-1,dat.n_cols);
    verbosity(Pa,"colsmooth: now do the smoothing of the row data by averaging every second row ",2,__FILE__,__LINE__);
    for (unsigned int rowo=0,rowi=0; rowo<out.n_rows; rowo+=2,rowi++)
    {
        //verbosity(Pa,"colsmooth: for rowi "+std::to_string(rowi),10,__FILE__,__LINE__);
        out.row(rowo)=dat.row(rowi);
        if(rowi>0 && rowi<dat.n_rows-1 && (rowo+1)<out.n_rows) out.row(rowo+1)=0.5*(dat.row(rowi-1)+dat.row(rowi+1)); //!< every second row is an interpolation of the neighbouring rows
    }

return out;

}

//%# Double data density by linear interpolation
inline arma::mat smooth(paramStruct Pa, arma::mat in)
{
    verbosity(Pa,"smooth: call colsmooth ",2,__FILE__,__LINE__);
    arma::mat out=colsmooth(in);
    verbosity(Pa,"smooth: smoothed column take transposed and smooth again ",2,__FILE__,__LINE__);
    out=colsmooth(out.t());
    verbosity(Pa,"smooth: return result ",2,__FILE__,__LINE__);
    return out;
}



#endif // SMOOTH_H
