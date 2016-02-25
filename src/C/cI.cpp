#include "cI.h"
#include "fft_minidft.h"

cI::cI(arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX)//,checkOperatorSize<T> *chkPointer)
{
std::cout << "setup L class" << std::endl;

this->_S=Sinp;
this->_outputM=arma::cx_mat(arma::prod(this->_S),Numwavfunc);

//this->_checkPointer=chkPointer;

myLatexClass=ltX;

std::string myFunctionString=" \\textbf{\\large{}forward fourier transformation operator cI} \\newline ";
myFunctionString+=" $I=\\mathbf{F}_{3}$ \\newline";
myFunctionString+=" uses fft to transform $\\Sigma_{k} S_{k} \\times N_{s}$ matrix from r to k space";
myFunctionString+=" ( forward == fftw backward transform)";

myLatexClass->setMyFunctionString(myFunctionString);

}

cI::~cI()
{
}

arma::cx_mat cI::operator*(arma::cx_mat  input)
{
arma::cx_mat outputMAT(input);
for(unsigned int i=0;i<input.n_cols;i++)
{
    //func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward)
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
    outputMAT.col(i)=func_fftn(input.col(i),this->_S,"backward");
}
return outputMAT;
}
arma::cx_mat cI::operator*(arma::cx_mat * input)
{
arma::cx_mat outputMAT(*input);
for(unsigned int i=0;i<input->n_cols;i++)
{
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
	outputMAT.col(i)=func_fftn(input->col(i),this->_S,"backward");
}
return outputMAT;
}




