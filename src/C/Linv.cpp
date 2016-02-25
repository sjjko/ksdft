#include "Linv.h"

Linv::Linv(const arma::mat G2inp,const arma::mat Rinp,const int Numwavfunc,latexComment *ltX)
//:checkOperatorSize<T>()
{
std::cout << "setup Linv class" << std::endl;
this->_G2i=G2inp;
this->_Ri=Rinp;
this->_Faktor=arma::mat(this->_G2i.n_rows,1,fill::zeros);
double detR=arma::det(Rinp);
for(int i=0;i<this->_Faktor.n_rows;i++)
{
    this->_Faktor(i,0)=-1./(G2inp(i,0)*detR);
}

this->_outputM=arma::cx_mat(G2inp.n_rows,Numwavfunc);
//verbosity("Linv constructor: set _Faktor 0 element to zero!")
this->_Faktor(0,0)=0;

myLatexClass=ltX;

std::string _myFunctionString=" \\textbf{\\large{}initiate inverse Laplacian operator class} \\newline ";
_myFunctionString+=" $L_{\\beta,\\alpha}^{-1}$ \\newline";
_myFunctionString+=" implemented for plane waves as: ";
_myFunctionString+=" $L^{-1} = {1 \\over det(\\mathbf{R})} D(G2)^{-1}$ ";
_myFunctionString+=" dim: $\\Sigma_{k} S_{k} \\times N_{s} $ ";

myLatexClass->setMyFunctionString(_myFunctionString);

}

Linv::~Linv()
{
}

/*applied to every column of input coefficient*/
arma::cx_mat  Linv::operator*(arma::cx_mat  input)
{
arma::cx_mat outputMAT(input);
for(unsigned int i=0;i<input.n_cols;i++)
{
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
	outputMAT.col(i)=this->_Faktor%input.col(i);
}
return outputMAT;
}
arma::cx_mat  Linv::operator*(arma::cx_mat * input)
{
arma::cx_mat outputMAT(*input);
for(unsigned int i=0;i<input->n_cols;i++)
{
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
	outputMAT.col(i)=this->_Faktor%(input->col(i));
}
return outputMAT;
}
/*arma::cx_vec  Linv::operator*(arma::cx_vec input)
{
	//check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
	return this->_Faktor%input;
}*/

//template class Linv<arma::mat>;
//template class Linv<arma::vec>;
//template class Linv<arma::cx_mat>;
