#ifndef OOP_H
#define OOP_H

//! the basevector product operator O
/*!
 This class provides the operator to compute the product of two basevectors
 in wavenumber space
 */

#include "main.h"
#include "latexComment.h"

//using namespace std;
//using namespace arma;
class latexComment;

class OOp
{
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream
        OOp(arma::mat Rinput,latexComment *ltX); //!< construct the overlap matrix class
        ~OOp();
        std::string getLatexString();
        inline void operator=(const OOp& other){this->_Ri=other._Ri;};  //!< assignment operator for overlap class
        inline OOp(const OOp& other){this->_Ri=other._Ri;}  //!< link overlap matrix classes
        arma::cx_mat  operator*(arma::cx_mat  input);  //!< operator (matrix) acts on (complex) state matrix
        arma::cx_mat  operator*(arma::cx_mat * input); //!< operator (matrix) acts on (complex) state matrix pointer
    private:
        arma::mat _Ri;  //!< the position of the cell as matrix of dimension $S_{k} \times 3$

};


#endif // O_H
