#include "cJ.h"
#include "fft_minidft.h"

cJ::cJ(arma::Col<double> Sinp,int Numwavfunc,latexComment *ltX)//,checkOperatorSize<T> *chkPointer)
{
std::cout << "setup L class" << std::endl;

this->_S=Sinp;
this->_outputM=arma::cx_mat(arma::prod(this->_S),Numwavfunc);

myLatexClass=ltX;

std::string myFunctionString=" \\textbf{\\large{}backward fourier transformation operator cJ} \\newline ";
myFunctionString+=" $I^{-1}=\\mathbf{F}_{3}^{-1}$ \\newline";
myFunctionString+=" uses fft to transform $\\Sigma_{k} S_{k} \\times N_{s}$ matrix from k to r space";
myFunctionString+=" ( backward == fftw forward transform)";

myLatexClass->setMyFunctionString(myFunctionString);
//this->_checkPointer=chkPointer;
}

cJ::~cJ()
{
}

arma::cx_mat cJ::operator*(arma::cx_mat  input)
{
arma::cx_mat outputMAT(input);
for(unsigned int i=0;i<input.n_cols;i++)
{
    //func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward)
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
    //this->_outputM.col(i)=1./arma::prod(this->_S)*func_fftn(input.col(i),this->_S,"forward");
    outputMAT.col(i)=1./arma::prod(this->_S)*func_fftn(input.col(i),this->_S,"forward");

}
return outputMAT;
//return this->_outputM;
}
arma::cx_mat cJ::operator*(arma::cx_mat * input)
{
arma::cx_mat outputMAT(*input);
for(unsigned int i=0;i<input->n_cols;i++)
{
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
	//this->_outputM.col(i)=1./arma::prod(this->_S)*func_fftn(input->col(i),this->_S,"forward");
	outputMAT.col(i)=1./arma::prod(this->_S)*func_fftn(input->col(i),this->_S,"forward");
}
return outputMAT;
//return this->_outputM;
}
arma::cx_mat cJ::operator*(arma::mat input)
{

//!< compute the FFT with real input matrix by copying it into the real part of a temporary complex matrix tmpM

arma::cx_mat tmpM(input,arma::zeros<arma::mat>(input.n_rows,input.n_cols));
arma::cx_mat outputMAT(tmpM);
for(unsigned int i=0;i<input.n_cols;i++)
{
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
	outputMAT.col(i)=1./arma::prod(this->_S)*func_fftn(tmpM.col(i),this->_S,"forward");
}
return outputMAT;
}

/*
arma::cx_vec cJ::operator*(arma::vec input)
{
arma::cx_vec cinput(input,arma::zeros<arma::vec>(input.n_elem));
                    //reinterpret_cast<arma::cx_vec> input;
return  1./arma::prod(this->_S)*func_fftn(cinput,this->_S,"forward");
}
*/

//template class L<arma::mat>;
//template class L<arma::vec>;
//template class L<arma::cx_mat>;



