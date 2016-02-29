#include "O.h"

OOp::OOp(arma::mat Rinput,latexComment *ltX)
{
    //! \brief initialize overlap matrix class
    //!
    //! \param Rinput: input matrix of grid point positions
    //! compute $O_{\beta,\alpha} = \int b_{\alpha}(\vec{x}) b_{\beta}(\vec{x}) dV$
    //! which for plane waves becomes:
    //! $O= det(\mathbf{R}) \mathbf{1}$


std::cout << "setup O class" << std::endl;
std::cout << "type:" << typeid(Rinput).name() << " Ri " << typeid(_Ri).name() << endl;
arma::mat Test(3,3);
_Ri=Test;
_Ri=Rinput;

myLatexClass=ltX;

std::string myFunctionString=" \\textbf{\\large{}initiate overlap matrix class} \\newline ";
myFunctionString+=" $O_{\\beta,\\alpha} = \\int b_{\\alpha}(\\vec{x}) b_{\\beta}(\\vec{x}) dV$ \\newline";
myFunctionString+=" implemented for plane waves: ";
myFunctionString+=" $O = det(\\mathbf{R}) \\mathbf{1}$";

myLatexClass->setMyFunctionString(myFunctionString);
//verbosity("Ri assigned !",2,__FILE__,__LINE__);
//verbosity("now go on !",2,__FILE__,__LINE__);
}


OOp::~OOp()
{
}

arma::cx_mat OOp::operator*(arma::cx_mat input)
{
return arma::det(_Ri)*input;
}
arma::cx_mat OOp::operator*(arma::cx_mat * input)
{
return arma::det(_Ri)* (*input);
}
/*
arma::cx_vec OOp::operator*(arma::cx_vec input)
{
return arma::det(_Ri)*input;
}
*/

//std::shared_ptr<arma::cx_mat> OOp::operator*(std::shared_ptr<arma::cx_mat> input)
//{
//bussardo
//return arma::det(this->Ri)* (*input.get());
//}

//template class O<arma::mat>;
//template class O<arma::vec>;
//template class O<arma::cx_mat>;
