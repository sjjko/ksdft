#include "L.h"

Lop::Lop(const arma::Mat<double> G2inp,const arma::Mat<double> G2CompInp,const arma::Mat<double> Rinp,const int Numwavfunc, latexComment *ltX)//,checkOperatorSize<T> *chkPointer)
{
std::cout << "setup L class" << std::endl;

this->_G2i=G2inp;
this->_Ri=Rinp;
this->_Faktor=-arma::det(Rinp)*G2inp;
this->_FaktorComp=-arma::det(Rinp)*G2CompInp;
_noActiveIndices=G2CompInp.n_rows;

this->_outputM=arma::cx_mat(G2inp.n_rows,Numwavfunc);

myLatexClass=ltX;

std::string _myFunctionString=" \\textbf{\\large{}initiate Laplacian operator class} \\newline ";
_myFunctionString+=" $L_{\\beta,\\alpha} = \\int b_{\\beta}(\\vec{x}) \\nabla^{2} b_{\\alpha}(\\vec{x}) dV$ \\newline";
_myFunctionString+=" implemented for plane waves as: ";
_myFunctionString+=" $L = -det(\\mathbf{R}) D(G2)$";
_myFunctionString+=" dim: $\\Sigma_{k} S_{k} \\times N_{s} $";

myLatexClass->setMyFunctionString(_myFunctionString);
//this->_checkPointer=chkPointer;
}

Lop::~Lop()
{
}

arma::cx_mat  Lop::operator*(arma::cx_mat  input)
{
arma::cx_mat outputMAT(input);
if(input.n_rows>_noActiveIndices){
    for(unsigned int i=0;i<input.n_cols;i++)
    {
        outputMAT.col(i)=this->_Faktor%input.col(i);
    }
}
else
{
    for(unsigned int i=0;i<input.n_cols;i++)
    {
        outputMAT.col(i)=this->_FaktorComp%input.col(i);
    }
}

return outputMAT;
}
arma::cx_mat  Lop::operator*(arma::cx_mat * input)
{
arma::cx_mat outputMAT(*input);
if(input->n_rows>_noActiveIndices){
    for(unsigned int i=0;i<input->n_cols;i++)
    {
        //check(_Faktor,input.col(i),"*"); // check if _Faktor has as many cols as input rows
        outputMAT.col(i)=this->_Faktor%(input->col(i));
    }
}
else
{
    for(unsigned int i=0;i<input->n_cols;i++)
    {
        outputMAT.col(i)=this->_FaktorComp%(input->col(i));
    }
}

return outputMAT;
}

//template class L<arma::mat>;
//template class L<arma::vec>;
//template class L<arma::cx_mat>;



