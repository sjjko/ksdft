#ifndef SLICE_H_INCLUDED
#define SLICE_H_INCLUDED

//using namespace std;
//using namespace arma;

#include "customAssert.h"

namespace myFunctions
{
inline arma::mat slice(arma::mat dat,arma::vec N,int n,int dir)
{

//! \brief
//!
//! \param N: dimensions of dat in a 3-vector
//! \param dat: 3d data set (any shape) of total size prod(N)=N(1)*N(2)*N(3)
//! \param dir: direction perpendicular to slice --- dir=0,1,2 gives yz,xz,yz planes
//! \param n: desired slice number from data; 1 <= n <= N(dir)
//! \return n-th dir-plane of dat (lower remaining dimension leading)
//!
//!

verbosity("slice: new out matrix ",2,__FILE__,__LINE__);

arma::mat out;

verbosity("start slicing ",2,__FILE__,__LINE__);

try
{

    //n=std::floor(n);
    if(n<1)
    {
        verbosity("slice does not exist ",2,__FILE__,__LINE__);
        //cout<<"Asking for non-existing slice!"<<endl;
        throw string("slice: Asking for non-existing slice!");
    }
    if(dir==2)
    {
        if(n>N[2])
        {
//cout << "Asking for non-existent slice" << n << endl;
            throw string("slice: Asking for non-existent slice!");
        }
//else reshape
        verbosity("data has dim "+std::to_string(N[0])+" , "+std::to_string(N[1])+" , "+std::to_string(N[2]) ,2,__FILE__,__LINE__);
        verbosity("slice: dir3 do reshape to "+std::to_string(N[0])+" x "+std::to_string(N[1])+" , "+std::to_string(N[2]) ,2,__FILE__,__LINE__);
        dat.reshape(N[0]*N[1],N[2]); //# Group into matrix with dir=3 as cols
        verbosity("take one column at n "+std::to_string(n) ,2,__FILE__,__LINE__);
        out=dat.col(n);
        verbosity("column has size "+std::to_string(out.n_elem) ,2,__FILE__,__LINE__);
        verbosity("slice: out do reshape to "+std::to_string(N[0])+" , "+std::to_string(N[1]) ,2,__FILE__,__LINE__);
        out.reshape(N[0],N[1]); // Take n-th col and reshape as slice
        verbosity("slice: dir3 finished reshaping" ,2,__FILE__,__LINE__);

    }
    else if(dir==1)
    {
        if(n>N[1])
        {
            //cout << "Asking for non-existent slice " << n << endl;
            throw string("slice: Asking for non-existent slice!");
        }
        verbosity("data has dim "+std::to_string(N[0])+" , "+std::to_string(N[1])+" , "+std::to_string(N[2]) ,2,__FILE__,__LINE__);
        verbosity("slice: dir2 do reshape to "+std::to_string(N[0])+" x "+std::to_string(N[1])+" , "+std::to_string(N[2]) ,2,__FILE__,__LINE__);
        dat.reshape(N[0]*N[1],N[2]); //%# Group to expose N[1]
        verbosity("slice: dir2 do an inplace transpose and convert to "+std::to_string(N[2])+" , "+std::to_string(N[1])+" x "+std::to_string(N[0]) ,2,__FILE__,__LINE__);
        inplace_trans(dat); //dat=//conj(dat'); //%# dat is now in order N[2],N[0]*N[1]
        verbosity("slice: dir2 reshape to "+std::to_string(N[2])+" x "+std::to_string(N[0])+" , "+std::to_string(N[1]) ,2,__FILE__,__LINE__);
        dat.reshape(N[2]*N[0],N[1]); //%# Form with dir=2 as cols
        verbosity("dir2: take one column at n "+std::to_string(n) ,2,__FILE__,__LINE__);
        out=dat.col(n);
        verbosity("column has size "+std::to_string(out.n_elem) ,2,__FILE__,__LINE__);
        verbosity("slice: out do reshape to "+std::to_string(N[2])+" , "+std::to_string(N[0]) ,2,__FILE__,__LINE__);
        out.reshape(N[2],N[0]); //%# Shape into slice
        verbosity("slice: dir2 do an inplace transpose and convert to "+std::to_string(N[0])+" , "+std::to_string(N[2]) ,2,__FILE__,__LINE__);
        inplace_trans(out);
        verbosity("slice: dir2 finished reshaping" ,2,__FILE__,__LINE__);
    }//out=conj(out'); //%# Reorder as N[0],N[2];
    else if(dir==0)
    {
        if(n>N[0])
        {
            throw string("slice: Asking for non-existent slice!");
        }
        verbosity("slice: dir1 do reshape to "+std::to_string(N[0])+" , "+std::to_string(N[1])+" x "+std::to_string(N[2]) ,2,__FILE__,__LINE__);
        dat.reshape(N[0],N[1]*N[2]); //%# Group to expose N[0]
        verbosity("slice: dir1 do an inplace transpose and convert to "+std::to_string(N[1])+" x "+std::to_string(N[2])+" , "+std::to_string(N[0]) ,2,__FILE__,__LINE__);
        inplace_trans(dat);//=conj(dat'); //%# dat is now N[1]*N[2],N[0]
        verbosity("dir1: take one column at n "+std::to_string(n) ,2,__FILE__,__LINE__);
        out=dat.col(n);
        verbosity("dir1 column has size "+std::to_string(out.n_elem) ,2,__FILE__,__LINE__);
        verbosity("slice: dir1 out do reshape to "+std::to_string(N[1])+" , "+std::to_string(N[2]) ,2,__FILE__,__LINE__);
        out.reshape(N[1],N[2]);
        verbosity("slice: dir1 finished reshaping" ,2,__FILE__,__LINE__);
    }
    else
    {
        throw string("Error in slice(): invalid choice for dir.  dir=");
//cout << "Error in slice(): invalid choice for dir.  dir=" << dir << endl;
    }
    verbosity("slice: out has dimensions "+std::to_string(out.size()) ,2,__FILE__,__LINE__);
    return out;
}
catch(string s)
{
    cout << s << endl;
}

return out;
}
}

#endif // SLICE_H_INCLUDED
