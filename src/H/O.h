#ifndef OOP_H
#define OOP_H

#include "main.h"
#include "latexComment.h"

//using namespace std;
//using namespace arma;

class OOp: public latexComment
{
    public:
        latexComment *myLatexClass; //!< public latex class, to be initialized in construction - used to give latex output stream
        OOp(arma::mat Rinput,latexComment *ltX); //!< construct the overlap matrix class
        ~OOp();
        std::string getLatexString();
       // OOp(const OOp& other); //copy constructor
        inline void operator=(const OOp& other){this->_Ri=other._Ri;};  //!< assignment operator for overlap class
        inline OOp(const OOp& other){this->_Ri=other._Ri;}  //!< link overlap matrix classes
        arma::cx_mat  operator*(arma::cx_mat  input);  //!< operator (matrix) acts on (complex) state matrix
        arma::cx_mat  operator*(arma::cx_mat * input); //!< operator (matrix) acts on (complex) state matrix pointer
       // arma::cx_vec  operator*(arma::cx_vec input);
        //std::shared_ptr<arma::cx_mat> operator*(std::shared_ptr<arma::cx_mat> input);
    private:
        arma::mat _Ri;  //!< the position of the cell as matrix of dimension $S_{k} \times 3$
        //arma::mat B(3,3);
//        std::auto_ptr<arma::mat> Rint;

};


#endif // O_H
