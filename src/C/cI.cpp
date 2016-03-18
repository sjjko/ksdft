#include "cI.h"
#include "fft_minidft.h"

using namespace myFunctions;

cI::cI( paramStruct Pa,arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX)//,checkOperatorSize<T> *chkPointer)
{
std::cout << "setup cI class" << std::endl;

this->_S=Sinp;
std::cout << "setup _outputM" << std::endl;

this->_outputM=arma::cx_mat(arma::prod(this->_S),Numwavfunc);
std::cout << "setup _Pa class" << std::endl;
this->_Pa=&Pa;

//this->_checkPointer=chkPointer;
std::cout << "setup Latex class" << std::endl;

myLatexClass=ltX;

std::string myFunctionString=" \\textbf{\\large{}forward fourier transformation operator cI} \\newline ";
myFunctionString+=" $I=\\mathbf{F}_{3}$ \\newline";
myFunctionString+=" uses fft to transform $\\Sigma_{k} S_{k} \\times N_{s}$ matrix from r to k space";
myFunctionString+=" ( forward == fftw backward transform)";

myLatexClass->setMyFunctionString(myFunctionString);
cout << "setup _uv matrix" << endl;
this->_uv=*this->_Pa->activeIndicesPtr;
cout << "setup full matrix" << endl;
this->_full=arma::zeros<arma::cx_mat>(arma::prod(_S),1);


}

cI::~cI()
{
}

arma::cx_mat cI::operator*(arma::cx_mat  input)
{


arma::cx_mat outputMAT(prod(_S),input.n_cols);


if(input.n_rows==prod(_S))
{
    for(unsigned int i=0;i<input.n_cols;i++)
    {
        //func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward)
        //check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
        outputMAT.col(i)=func_fftn(input.col(i),this->_S,"backward");
    }
}
else
{
   // verbosity(*_Pa,"CI: do the transform",2,__FILE__,__LINE__);

    for(unsigned int i=0;i<input.n_cols;i++)
    {
    //verbosity(*_Pa,"CI: do it for column"+std::to_string(i),2,__FILE__,__LINE__);

        //func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward)
        //check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
//        cout << size(full(_Pa->activeIndices)) << endl;
//                cout << size(input) << endl;
//                cout << size(input.col(i)) << endl;
//                cout << size(test) << endl;

//        for(int j=0;j<_Pa->numberOfActiveIndices;j++)
//        {
//        cout << uv(j) << " " << j << endl;
//        full(uv(j),0)=input(j,i);
//        }
        _full=arma::zeros<arma::cx_mat>(arma::prod(_S),1);
        _full(_uv)=input.col(i);
        outputMAT.col(i)=func_fftn(_full,this->_S,"backward");
    }
}

return outputMAT;
}
arma::cx_mat cI::operator*(arma::cx_mat * input)
{
//arma::cx_mat outputMAT(*input);
arma::cx_mat outputMAT(prod(_S),input->n_cols);


    //verbosity(*_Pa,"CI: do the * transform",2,__FILE__,__LINE__);

if(input->n_rows==prod(_S))
{
    for(unsigned int i=0;i<input->n_cols;i++)
    {
        //check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
        outputMAT.col(i)=func_fftn(input->col(i),this->_S,"backward");
    }
    }
else
{
    for(unsigned int i=0;i<input->n_cols;i++)
    {
        //verbosity(*_Pa,"CI: do * for column"+std::to_string(i),2,__FILE__,__LINE__);

        //func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward)
        //check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
        _full=arma::zeros<arma::cx_mat>(arma::prod(_S),1);
        _full(_uv)=input->col(i);
        outputMAT.col(i)=func_fftn(_full,this->_S,"backward");
    }
}

return outputMAT;
}




