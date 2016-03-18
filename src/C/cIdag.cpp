#include "cIdag.h"
#include "fft_minidft.h"

cIdag::cIdag(const paramStruct Pa,arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX)//,checkOperatorSize<T> *chkPointer)
{
std::cout << "setup cIdag class" << std::endl;

this->_S=Sinp;
this->_outputM=arma::cx_mat(arma::prod(this->_S),Numwavfunc);
this->_Pa=&Pa;

myLatexClass=ltX;

std::string  myFunctionString=" \\textbf{\\large{}hermitian conjugate forward fourier transformation operator cIdag} \\newline ";
myFunctionString+=" $I^{\\dagger}=\\mathbf{F}^{\\dagger}_{3}$ \\newline";
myFunctionString+=" uses fft to transform $\\Sigma_{k} S_{k} \\times N_{s}$ matrix from r to k space";
myFunctionString+=" (conjugated forward == fftw forward transform)";

myLatexClass->setMyFunctionString(myFunctionString);

this->_uv=*this->_Pa->activeIndicesPtr;

//this->_checkPointer=chkPointer;
}

cIdag::~cIdag()
{
}

arma::cx_mat cIdag::operator*(arma::cx_mat  input)
{
//arma::cx_mat outputMAT(input);

arma::cx_mat outputMAT=zeros<cx_mat>(_Pa->numberOfActiveIndices,input.n_cols);
arma::cx_mat full(input.n_rows,1);

for(unsigned int i=0;i<input.n_cols;i++)
{
    //func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward)
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
    //this->_outputM.col(i)=func_fftn(input.col(i),this->_S,"forward");
    full=func_fftn(input.col(i),this->_S,"forward");
    outputMAT.col(i)=full(_uv);

}
return outputMAT;
}
arma::cx_mat cIdag::operator*(arma::cx_mat * input)
{
//arma::cx_mat outputMAT(*input);

arma::cx_mat outputMAT=zeros<cx_mat>(_Pa->numberOfActiveIndices,input->n_cols);
arma::cx_mat full(input->n_rows,1);

for(unsigned int i=0;i<input->n_cols;i++)
{
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
	full=func_fftn(input->col(i),this->_S,"forward");
	outputMAT.col(i)=full(_uv);
}
return outputMAT;
}

//template class L<arma::mat>;
//template class L<arma::vec>;
//template class L<arma::cx_mat>;



