#include "cJdag.h"
#include "fft_minidft.h"

cJdag::cJdag(arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX)//,checkOperatorSize<T> *chkPointer)
{
std::cout << "setup L class" << std::endl;

this->_S=Sinp;
this->_outputM=arma::cx_mat(arma::prod(this->_S),Numwavfunc);

myLatexClass=ltX;

std::string myFunctionString=" \\textbf{\\large{}hermitian conjugate backward fourier transformation operator cIdag} \\newline ";
myFunctionString+=" $(I^{\\dagger})^{-1}=\\mathbf{F}^{\\dagger}_{3}$ \\newline";
myFunctionString+=" uses fft to transform $\\Sigma_{k} S_{k} \\times N_{s}$ matrix from k to r space";
myFunctionString+=" (conjugated forward == fftw forward transform)";

myLatexClass->setMyFunctionString(myFunctionString);
//this->_checkPointer=chkPointer;
}

cJdag::~cJdag()
{
}

arma::cx_mat cJdag::operator*(arma::cx_mat  input)
{
arma::cx_mat outputMAT(input);
for(unsigned int i=0;i<input.n_cols;i++)
{
    //func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward)
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
    //this->_outputM.col(i)=1./arma::prod(this->_S)*func_fftn(input.col(i),this->_S,"backward");
    outputMAT.col(i)=1./arma::prod(this->_S)*func_fftn(input.col(i),this->_S,"backward");

}
return outputMAT;
}
arma::cx_mat cJdag::operator*(arma::cx_mat * input)
{
arma::cx_mat outputMAT(*input);
for(unsigned int i=0;i<input->n_cols;i++)
{
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
	this->_outputM.col(i)=1./arma::prod(this->_S)*func_fftn(input->col(i),this->_S,"backward");
}
return outputMAT;
}

//template class L<arma::mat>;
//template class L<arma::vec>;
//template class L<arma::cx_mat>;



